#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Helper: Kernel Function (Anisotropic)
arma::mat cov_kernel_gau(arma::mat x1, arma::mat x2, arma::vec theta) {
  arma::uword n1 = x1.n_rows; arma::uword n2 = x2.n_rows; arma::uword d = x1.n_cols;
  arma::mat K(n1, n2);
  for(arma::uword i=0; i<n1; ++i) {
    for(arma::uword j=0; j<n2; ++j) {
      double w_dist = 0;
      for(arma::uword k=0; k<d; ++k) {
        double diff = x1(i,k) - x2(j,k);
        w_dist += (diff * diff) / theta(k);
      }
      K(i,j) = exp(-w_dist);
    }
  }
  return K;
}

arma::mat cov_kernel_m32(arma::mat x1, arma::mat x2, arma::vec theta) {
  arma::uword n1 = x1.n_rows; arma::uword n2 = x2.n_rows; arma::uword d_cols = x1.n_cols;
  arma::mat K(n1, n2);
  double sqrt3 = std::sqrt(3.0);

  for(arma::uword i=0; i<n1; ++i) {
    for(arma::uword j=0; j<n2; ++j) {
      double prod = 1.0; // Initialize product
      for(arma::uword k=0; k<d_cols; ++k) {
        double d_k = std::abs(x1(i,k) - x2(j,k)) / theta(k);
        // Apply Matern formula to each dimension k
        prod *= (1.0 + sqrt3 * d_k) * std::exp(-sqrt3 * d_k);
      }
      K(i,j) = prod;
    }
  }
  return K;
}

arma::mat cov_kernel_m52(arma::mat x1, arma::mat x2, arma::vec theta) {
  arma::uword n1 = x1.n_rows; arma::uword n2 = x2.n_rows; arma::uword d_cols = x1.n_cols;
  arma::mat K(n1, n2);
  double sqrt5 = std::sqrt(5.0);

  for(arma::uword i=0; i<n1; ++i) {
    for(arma::uword j=0; j<n2; ++j) {
      double prod = 1.0;
      for(arma::uword k=0; k<d_cols; ++k) {
        double d_k = std::abs(x1(i,k) - x2(j,k)) / theta(k);
        // Apply Matern formula to each dimension k
        prod *= (1.0 + sqrt5 * d_k + (5.0/3.0) * d_k * d_k) * std::exp(-sqrt5 * d_k);
      }
      K(i,j) = prod;
    }
  }
  return K;
}

// [[Rcpp::export]]
arma::mat compute_kernel(arma::mat x1, arma::mat x2, arma::vec theta, std::string type) {
  if (type == "Gaussian") {
    return cov_kernel_gau(x1, x2, theta);
  } else if (type == "Matern3_2") {
    return cov_kernel_m32(x1, x2, theta);
  } else if (type == "Matern5_2") {
    return cov_kernel_m52(x1, x2, theta);
  } else {
    Rcpp::stop("Unknown kernel type");
    return mat();
  }
}

// [[Rcpp::export]]
Rcpp::List gp_predict_cpp(arma::mat X_new, arma::mat X, arma::vec Y,
                          std::string Cov_fun, Rcpp::List Parameters) {

  arma::vec theta = Parameters["theta"];
  double g = Parameters["g"];
  double sigma2 = Parameters["sigma2"];
  double beta0 = Parameters["beta0"];

  // 1. Build and solve the training system
  arma::mat K = compute_kernel(X, X, theta, Cov_fun);
  K.diag() += g;

  // 2. Compute the cross-kernel (transpose it once for better memory access)
  arma::mat K_star_T = compute_kernel(X_new, X, theta, Cov_fun).t();

  // 3. Solve systems: K * alpha = (Y - beta0) and K * B = K_star_T
  arma::vec alpha = arma::solve(K, Y - beta0, arma::solve_opts::fast);
  arma::mat B = arma::solve(K, K_star_T, arma::solve_opts::fast);

  // 4. Mean prediction: mu = beta0 + K_star * alpha
  arma::vec mu = beta0 + (K_star_T.t() * alpha);

  // 5. Variance prediction
  arma::vec s2(X_new.n_rows);
  for(arma::uword i = 0; i < X_new.n_rows; ++i) {
    s2(i) = sigma2 * ((1.0 + g) - arma::dot(K_star_T.col(i), B.col(i)));
  }

  return Rcpp::List::create(
    Rcpp::Named("mean") = mu,
    Rcpp::Named("sd2") = s2
  );
}

arma::mat inverse_update_cpp(arma::mat Rs_inv, arma::vec rs_add, double g) {
  arma::uword d = Rs_inv.n_cols;

  // A = Rs_inv * rs_add
  arma::vec A = Rs_inv * rs_add;

  // Q = (1 + g) - rs_add' * A
  double Q = (1.0 + g) - as_scalar(rs_add.t() * A);

  arma::mat Rs_inv_next(d + 1, d + 1);
  arma::vec temp = A / Q;

  // Fill the blocks
  Rs_inv_next.submat(0, 0, d - 1, d - 1) = Rs_inv + (temp * A.t());
  Rs_inv_next.submat(0, d, d - 1, d)     = -temp;
  Rs_inv_next.submat(d, 0, d, d - 1)     = -temp.t();
  Rs_inv_next.at(d, d)                   = 1.0 / Q;

  return Rs_inv_next;
}


// [[Rcpp::export]]
Rcpp::List GEMSS_cpp_update_sig(arma::mat X, arma::vec Y, arma::mat X_val, arma::vec Y_val, arma::uvec initial_idx, arma::uword ns, double c1, arma::uword n_srs, arma::uword n_top,
                     arma::vec theta, double nugget, double beta0, std::string Cov_fun, bool print_result) {

  arma::uword N = X.n_rows;
  arma::uvec index(ns);
  arma::vec r_sq_history(ns);
  double sig2;

  arma::uword current_size = initial_idx.n_elem;
  index.head(current_size) = initial_idx - 1; // Fill the start

  // PRE-ALLOCATE instead of join_cols
  arma::mat X_sub_all(ns, X.n_cols);
  arma::vec Y_sub_all(ns);
  X_sub_all.rows(0, current_size - 1) = X.rows(index.head(current_size));
  Y_sub_all.head(current_size) = Y.rows(index.head(current_size));

  // --- THE FAST MASKING START ---
  std::vector<bool> is_used(N, false);

  // Only mark the ACTUAL initial points
  for(arma::uword i = 0; i < current_size; ++i) {
    arma::uword idx = index(i);
    if(idx < N) is_used[idx] = true;
  }

  // Build the initial screening set efficiently
  std::vector<uword> initial_scr_vec;
  initial_scr_vec.reserve(N - current_size);
  for(arma::uword i = 0; i < N; ++i) {
    if(!is_used[i]) {
      initial_scr_vec.push_back(i);
    }
  }
  arma::uvec scr_ind = conv_to<uvec>::from(initial_scr_vec);


  // Initial full inverse (only done once)
  arma::mat X_sub = X_sub_all.rows(0, current_size - 1);
  arma::vec Y_sub_c = Y_sub_all.head(current_size) - beta0;
  arma::mat K_ff = compute_kernel(X_sub, X_sub, theta, Cov_fun);
  K_ff.diag() += nugget;
  arma::mat K_inv = inv(K_ff);

  // alpha
  arma::vec alpha = K_inv * Y_sub_c;
  sig2 = arma::as_scalar(Y_sub_c.t() * alpha) / (double)(current_size);

  arma::mat K_val_f(X_val.n_rows, ns);
  K_val_f.cols(0, current_size - 1) = compute_kernel(X_val, X_sub, theta, Cov_fun);
  arma::vec mu_val = (K_val_f.cols(0, current_size-1) * alpha) + beta0;
  double mst = mean(square(Y_val - mean(Y)));
  // double mst = mean(square(Y_val - mu_val));

  while(current_size < ns) {
    // 1. PREDICTION & VARIANCE (O(n^2))
    arma::mat X_scr = X.rows(scr_ind);
    arma::mat K_sf = compute_kernel(X_scr, X_sub, theta, Cov_fun);
    arma::vec mu_s = (K_sf * alpha) + beta0;
    arma::vec var_s(scr_ind.n_elem);
    for(arma::uword i=0; i < scr_ind.n_elem; ++i) {
      arma::rowvec k_i = K_sf.row(i);
      var_s(i) = sig2 * ((1.0 + nugget) - as_scalar(k_i * K_inv * k_i.t()));
    }

    // 2. SELECTION
    arma::vec d_y = square(Y.rows(scr_ind) - mu_s);
    arma::vec Gu;
    if (c1 == -1.0) {
      arma::vec c1_vec = 1.0 / (2.0 * (square(d_y) / (3.0 * square(var_s))) + 1.0);
      Gu = c1_vec % d_y + (1.0 - c1_vec) % var_s;
    } else {
      Gu = c1 * d_y + (1.0 - c1) * var_s;
    }
    arma::uvec sorted_local_idx = sort_index(Gu, "descend");
    arma::uword best_global_idx = scr_ind(sorted_local_idx(0));

    index(current_size) = best_global_idx;

    // 3. INCREMENTAL UPDATE (Update model with the new point)
    arma::mat k_new_mat = compute_kernel(X_sub, X.row(best_global_idx), theta, Cov_fun);
    K_inv = inverse_update_cpp(K_inv, k_new_mat.col(0), nugget);


    X_sub_all.row(current_size) = X.row(best_global_idx);
    Y_sub_all(current_size) = Y(best_global_idx);
    X_sub = X_sub_all.rows(0, current_size);
    Y_sub_c = Y_sub_all.head(current_size+1) - beta0;

    // 4. PREDICTIVE R-SQUARE (Evaluate the NEWLY updated model)
    // Recalculate alpha for the expanded system
    alpha = K_inv * Y_sub_c;
    sig2 = arma::as_scalar(Y_sub_c.t() * alpha) / (double)(current_size + 1);

    K_val_f.col(current_size) = compute_kernel(X_val, X.row(best_global_idx), theta, Cov_fun);
    mu_val = (K_val_f.cols(0, current_size) * alpha) + beta0;
    double mse = mean(square(Y_val - mu_val));
    r_sq_history(current_size) = 1.0 - (mse / mst); // Save at the current size index

    // 5. REGENERATE SCREENED CANDIDATES
    arma::uword actual_top = std::min((uword)n_top, (uword)sorted_local_idx.n_elem - 1);
    arma::uvec top_cand = scr_ind.elem(sorted_local_idx.subvec(1, actual_top));
    arma::uvec srs_ind = randi<uvec>(n_srs, distr_param(0, N - 1));
    arma::uvec combined = unique(join_cols(top_cand, srs_ind));

    // Create a boolean mask of all points currently in the active set
    std::vector<bool> in_active_set(N, false);
    for(arma::uword i = 0; i < current_size + 1; ++i) {
      in_active_set[index(i)] = true;
    }

    std::vector<uword> next_scr_vec;
    next_scr_vec.reserve(combined.n_elem); // Pre-allocate for speed

    for(arma::uword i = 0; i < combined.n_elem; ++i) {
      // Instant look-up instead of searching through 'index'
      if(!in_active_set[combined(i)]) {
        next_scr_vec.push_back(combined(i));
      }
    }
    scr_ind = conv_to<uvec>::from(next_scr_vec);


    if (print_result) {
      Rcpp::Rcout << "Size: " << current_size + 1
                  << " | R2: " << std::fixed << std::setprecision(4) << r_sq_history(current_size)
                  << std::endl;
    }
    current_size++;
  }

  return Rcpp::List::create(
    Rcpp::_["index"] = index + 1, // Back to 1-based for R
    Rcpp::_["r_sq"] = r_sq_history,
    Rcpp::_["sigma2"] = sig2
  );
}

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

arma::mat inverse_update_cpp(arma::mat U, const arma::vec& rs_add, double g) {
  arma::uword d = U.n_cols;
  double tol = 1e-12;

  // Solve t(U) * v = rs_add
  arma::vec v = arma::solve(arma::trimatu(U).t(), rs_add, arma::solve_opts::fast);

  // Schur complement (alpha^2)
  double alpha2 = 1.0 + g - arma::dot(v, v);

  if (alpha2 < -tol) {
    Rcpp::stop("Cholesky update failed: Matrix is not numerically SPD.");
  }

  double alpha = std::sqrt(std::max(alpha2, 1e-14)); // Jitter for stability

  // Build the new matrix explicitly
  arma::mat U_next(d + 1, d + 1, arma::fill::zeros);
  U_next.submat(0, 0, d - 1, d - 1) = U;
  U_next.submat(0, d, d - 1, d) = v;
  U_next(d, d) = alpha;

  return U_next;
}


// [[Rcpp::export]]
Rcpp::List GEMSS_cpp_update_sig(const arma::mat& X, const arma::vec& Y, const arma::mat& X_val, const arma::vec& Y_val, const arma::uvec& initial_idx,
                                arma::uword ns, double c1, arma::uword n_srs, arma::uword n_top, const arma::vec& theta, double nugget, double beta0,
                                const std::string& Cov_fun, bool print_result) {
  arma::uword N = X.n_rows;
  arma::uword p = X.n_cols;

  arma::uvec index(ns);
  arma::vec r_sq_history(ns, arma::fill::zeros);

  arma::uword current_size = initial_idx.n_elem;
  index.head(current_size) = initial_idx - 1;

  arma::mat X_sub_all(ns, p);
  arma::vec Y_sub_c(ns);

  X_sub_all.rows(0, current_size - 1) = X.rows(index.head(current_size));
  Y_sub_c.head(current_size)   = Y.elem(index.head(current_size)) - beta0;

  // persistent active mask
  std::vector<unsigned char> is_active(N, 0);
  for (arma::uword i = 0; i < current_size; ++i) {
    is_active[index(i)] = 1;
  }

  std::vector<arma::uword> scr_vec;
  scr_vec.reserve(N - current_size);
  for (arma::uword i = 0; i < N; ++i) {
    if (!is_active[i]) scr_vec.push_back(i);
  }
  arma::uvec scr_ind = arma::conv_to<arma::uvec>::from(scr_vec);

  // initial model
  arma::mat K_ff = compute_kernel(
    X_sub_all.rows(0, current_size - 1),
    X_sub_all.rows(0, current_size - 1),
    theta, Cov_fun
  );
  K_ff.diag() += nugget;
  arma::mat U = arma::chol(K_ff);

  arma::vec w, alpha;
  w = arma::solve(arma::trimatu(U).t(), Y_sub_c.head(current_size));
  alpha = arma::solve(arma::trimatu(U), w);
  double sig2 = arma::as_scalar(Y_sub_c.head(current_size).t() * alpha) / double(current_size);

  arma::mat K_val_f(X_val.n_rows, ns, arma::fill::none);
  K_val_f.cols(0, current_size - 1) =
    compute_kernel(X_val, X_sub_all.rows(0, current_size - 1), theta, Cov_fun);

  arma::vec mu_val(X_val.n_rows, arma::fill::none);

  double mst = arma::mean(arma::square(Y_val - arma::mean(Y)));
  double mse = 0.0;

  // reusable workspace
  arma::mat K_sf, V;
  arma::vec mu_s, v_norms, var_s, y_scr, d_y, Gu;
  std::vector<arma::uword> next_scr_vec;
  next_scr_vec.reserve(n_top + n_srs);

  std::vector<unsigned char> in_next_scr(N, 0);

  while (current_size < ns) {
    K_sf = compute_kernel(X.rows(scr_ind), X_sub_all.rows(0, current_size - 1), theta, Cov_fun);

    mu_s = K_sf * alpha;
    mu_s += beta0;

    V = arma::solve(arma::trimatu(U).t(), K_sf.t(), arma::solve_opts::fast);
    v_norms = arma::sum(arma::square(V), 0).t();

    var_s = (1.0 + nugget) - v_norms;
    var_s *= sig2;

    y_scr = Y.elem(scr_ind);
    d_y = y_scr - mu_s;
    d_y %= d_y;

    if (c1 == -1.0) {
      arma::vec tmp1 = d_y;
      tmp1 %= tmp1;

      arma::vec tmp2 = var_s;
      tmp2 %= tmp2;
      tmp2 *= 3.0;

      tmp1 /= tmp2;
      tmp1 *= 2.0;
      tmp1 += 1.0;

      arma::vec c1_vec = 1.0 / tmp1;
      Gu = c1_vec % d_y + (1.0 - c1_vec) % var_s;
    } else {
      Gu = c1 * d_y + (1.0 - c1) * var_s;
    }

    arma::uvec sorted_local_idx = arma::sort_index(Gu, "descend");
    arma::uword best_global_idx = scr_ind(sorted_local_idx(0));

    index(current_size) = best_global_idx;
    is_active[best_global_idx] = 1;

    X_sub_all.row(current_size) = X.row(best_global_idx);
    Y_sub_c(current_size) =  Y(best_global_idx) - beta0;

    arma::vec k_new = compute_kernel(
      X_sub_all.rows(0, current_size - 1),
      X.row(best_global_idx),
      theta, Cov_fun
    ).col(0);

    if ((current_size + 1) % 50 == 0) {
      arma::mat K_ff_direct = compute_kernel(
        X_sub_all.rows(0, current_size),
        X_sub_all.rows(0, current_size),
        theta, Cov_fun
      );
      K_ff_direct.diag() += nugget;
      U = arma::chol(K_ff_direct);
    } else {
      U = inverse_update_cpp(U, k_new, nugget);
    }

    w = arma::solve(arma::trimatu(U).t(), Y_sub_c.head(current_size + 1));
    alpha = arma::solve(arma::trimatu(U), w);
    sig2 = arma::as_scalar(Y_sub_c.head(current_size + 1).t() * alpha) / double(current_size + 1);

    K_val_f.col(current_size) = compute_kernel(X_val, X.row(best_global_idx), theta, Cov_fun);
    mu_val = K_val_f.cols(0, current_size) * alpha;
    mu_val += beta0;

    mse = arma::mean(arma::square(Y_val - mu_val));
    r_sq_history(current_size) = 1.0 - mse / mst;

    next_scr_vec.clear();
    next_scr_vec.reserve(n_top + n_srs);

    if (sorted_local_idx.n_elem > 1) {
      arma::uword actual_top = std::min<arma::uword>(n_top, sorted_local_idx.n_elem - 1);
      for (arma::uword j = 1; j <= actual_top; ++j) {
        arma::uword idx = scr_ind(sorted_local_idx(j));
        if (!is_active[idx] && !in_next_scr[idx]) {
          next_scr_vec.push_back(idx);
          in_next_scr[idx] = 1;
        }
      }
    }

    for (arma::uword j = 0; j < n_srs; ++j) {
      arma::uword idx = arma::randi<arma::uword>(arma::distr_param(0, N - 1));
      if (!is_active[idx] && !in_next_scr[idx]) {
        next_scr_vec.push_back(idx);
        in_next_scr[idx] = 1;
      }
    }

    scr_ind = arma::conv_to<arma::uvec>::from(next_scr_vec);

    for (arma::uword idx : next_scr_vec) {
      in_next_scr[idx] = 0;
    }

    if (print_result) {
      Rcpp::Rcout << "Size: " << current_size + 1
                  << " | R2: " << std::fixed << std::setprecision(4)
                  << r_sq_history(current_size) << "\n";
    }

    ++current_size;
  }

  return Rcpp::List::create(
    Rcpp::_["index"] = index + 1,
    Rcpp::_["r_sq"] = r_sq_history,
    Rcpp::_["sigma2"] = sig2
  );
}

// [[Rcpp::export]]
Rcpp::List gemss_removal_cpp(arma::mat X, arma::vec Y, const arma::vec& theta, double nugget, double beta0, double sig2,
                             std::string Cov_fun, int n_remove, double c1) {
  int n = X.n_rows;
  arma::uvec removed_order(n_remove);
  arma::uvec current_indices = arma::regspace<arma::uvec>(0, n - 1);

  arma::mat K = compute_kernel(X, X, theta, Cov_fun);
  K.diag() += nugget;
  arma::mat R_inv = arma::inv_sympd(K);

  arma::vec Y_centered = Y - beta0;
  arma::vec inv_diag, eps_y, eps_r, cri, b;

  for (int j = 0; j < n_remove; ++j) {
    arma::uword n_cur = current_indices.n_elem;

    // Dubrule Formula
    inv_diag = R_inv.diag();
    eps_y = (R_inv * Y_centered) / inv_diag;
    eps_r = sig2 / inv_diag;

    // Compute GEMSS Criterion for each point
    cri.set_size(n_cur);
    if (c1 == -1.0) {
      for (arma::uword i = 0; i < n_cur; ++i) {
        double var_r = eps_r(i);
        double err_y = eps_y(i);

        double err_y2 = err_y * err_y;
        double err_y4 = err_y2 * err_y2;
        double var_r2 = var_r * var_r;

        double c1 = 1.0 / (2.0 * (err_y4 / (3.0 * var_r2)) + 1.0);
        cri(i) = c1 * err_y2 + (1.0 - c1) * var_r;
      }
    } else {
      for (arma::uword i = 0; i < n_cur; ++i) {
        double var_r = eps_r(i);
        double err_y = eps_y(i);

        cri(i) = c1 * err_y * err_y + (1.0 - c1) * var_r;
      }

    }


    arma::uword loc = cri.index_min();
    removed_order(j) = current_indices(loc) + 1;

    double d_val = inv_diag(loc);

    // update
    // O(n^2) In-place Inverse Update
    b = R_inv.col(loc);
    R_inv -= (b * b.t()) / d_val;

    R_inv.shed_col(loc);
    R_inv.shed_row(loc);

    Y_centered.shed_row(loc);
    current_indices.shed_row(loc);
  }

  return Rcpp::List::create(
    Rcpp::Named("removed_order") = removed_order,
    Rcpp::Named("index") = current_indices + 1
  );
}


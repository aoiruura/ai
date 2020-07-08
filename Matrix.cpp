#include <iostream>
#include <vector>
#include <queue>
#include <cmath>
#include <iomanip>
#include <algorithm>
using namespace std;

using D = double;
const D EPS = 1e-10;

template<class T> struct Matrix {
	vector<vector<T> > val;
	Matrix(int n, int m, T x = 0) : val(n, vector<T>(m, x)) {}
	void init(int n, int m, T x = 0) { val.assign(n, vector<T>(m, x)); }
	size_t size() const { return val.size(); }
	inline vector<T>& operator [](int i) { return val[i]; }
};

template<class T> ostream& operator << (ostream& c, Matrix<T> A) {
	c << endl;
	for (int i = 0; i < A.size(); ++i) {
		for (int j = 0; j < A[i].size(); ++j) {
			c << A[i][j] << (j == A[i].size() - 1 ? "\n" : " ");
		}
	}
	return c;
}

template<class T> istream& operator >> (istream& c, Matrix<T>& A) {
	for (int i = 0; i < A.size(); ++i) {
		for (int j = 0; j < A[0].size(); ++j) {
			c >> A[i][j];
		}
	}
	return c;
}

template<class T> Matrix<T> operator * (Matrix<T> A, Matrix<T> B) {
	Matrix<T> ret(A.size(), B[0].size());
	for (int i = 0; i < A.size(); ++i) {
		for (int j = 0; j < B[0].size(); ++j) {
			for (int k = 0; k < B.size(); ++k) {
				ret[i][j] = A[i][k] * B[k][j];
			}
		}
	}
	return ret;
}

template<class T> Matrix<T> operator - (Matrix<T> A, Matrix<T> B) {
	Matrix<T> ret(A.size(), A[0].size());
	for (int i = 0; i < A.size(); ++i) {
		for (int j = 0; j < A[0].size(); ++j) {
			ret[i][j] = A[i][j] - B[i][j];
		}
	}
	return ret;
}

template<class T> Matrix<T> operator + (Matrix<T> A, Matrix<T> B) {
	Matrix<T> ret(A.size(), A[0].size());
	for (int i = 0; i < A.size(); ++i) {
		for (int j = 0; j < A[0].size(); ++j) {
			ret[i][j] = A[i][j] + B[i][j];
		}
	}
	return ret;
}

template<class T> void operator * (T t, Matrix<T>& A) {
	for (int i = 0; i < A.size(); ++i) {
		for (int j = 0; j < A[0].size(); ++j) {
			A[i][j] *= t;
		}
	}
}

template<class T> int GaussJordan(Matrix<T>& A, bool is_extended = false) {
	int m = A.size(), n = A[0].size();
	int rank = 0;
	for (int col = 0; col < n; ++col) {
		if (is_extended && col == n - 1) break;

		int pivot = -1;
		T ma = EPS;
		for (int row = rank; row < m; ++row) {
			if (abs(A[row][col]) > ma) {
				ma = abs(A[row][col]);
				pivot = row;
			}
		}

		if (pivot == -1) continue;

		swap(A[pivot], A[rank]);

		auto fac = A[rank][col];
		for (int col2 = 0; col2 < n; ++col2) A[rank][col2] /= fac;

		for (int row = 0; row < m; ++row) {
			if (row != rank && abs(A[row][col]) > EPS) {
				auto fac = A[row][col];
				for (int col2 = 0; col2 < n; ++col2) {
					A[row][col2] -= A[rank][col2] * fac;
				}
			}
		}
		++rank;
	}
	return rank;
}

template<class T> vector<T> liner_equation(Matrix<T> A, vector<T> b) {
	int m = A.size(), n = A[0].size();
	Matrix<T> M(m, n + 1);
	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < n; ++j) M[i][j] = A[i][j];
		M[i][n] = b[i];
	}
	int rank = GaussJordan(M, true);

	vector<T> res;
	for (int row = 0; row < m; ++row) if (abs(M[row][n]) > EPS) return res;

	res.assign(n, 0);
	for (int i = 0; i < rank; ++i) res[i] = M[i][n];
	return res;
}

template<class T> Matrix<T> Judge_regular(Matrix<T>& A) {
	int m = A.size(), n = A[0].size();
	int delta = 0;
	//[A, E]の形を作る
	for (int row = 0; row < m; ++row) {
		for (int col = 0; col < n; ++col) {
			if (col == delta) A[row].push_back(1);
			else A[row].push_back(0);
		}
		delta++;
	}
	//行は2倍の長さになるがさがすのはnまで
	//n *= 2;
	int rank = 0;
	for (int col = 0; col < n; ++col) {
		
		int pivot = -1;
		T ma = EPS;
		for (int row = rank; row < m; ++row) {
			if (abs(A[row][col]) > ma) {
				ma = abs(A[row][col]);
				pivot = row;
			}
		}
		if (pivot == -1) continue;
		swap(A[pivot], A[rank]);

		auto fac = A[rank][col];
		for (int col2 = 0; col2 < 2 * n; ++col2) A[rank][col2] /= fac;
		for (int row = 0; row < m; ++row) {
			if (row != rank && abs(A[row][col]) > EPS) {
				auto fac = A[row][col];
				for (int col2 = 0; col2 < 2 * n; ++col2) {
					A[row][col2] -= A[rank][col2] * fac;
				}
			}
		}
		++rank;
	}
	Matrix<T> ret(m, n, 0);
	if (rank == m) {
		for (int row = 0; row < m; ++row) {
			for (int col = n; col < 2 * n; ++col) {
				ret[row][col - n] = A[row][col];
			}
		}
		return ret;
	}
	else return ret;
}

int main() {
	int M, N;
	cin >> M >> N;
	Matrix<D> Mat(M, N);
	cin >> Mat;
	auto TMat = Judge_regular(Mat);
	cout << TMat;
};

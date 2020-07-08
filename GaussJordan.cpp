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

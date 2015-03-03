#include <iostream>
#include <vector>
#include <set>
#include <cmath>
#include <cstdlib>
using namespace std;
#define EPS 1E-9
#define DEBUG 0
int count;
inline int identity_col (const vector <vector <double> > & A, int c) {
  int count = 0, row;
  for (int r=0; r<A.size(); r++) 
    if (A[r][c] > EPS) { count++; row = r; }  
  return (count == 1) ? row : -1;
}
void canonicalize ( vector <vector <double> > & A,
		    vector <double>& B, 
		    vector <double>& C,
		    vector <int>& BasicVarR,   // basic variable of each row
		    double & obj                // objective value
		    )
{
  int m = A.size(), n = C.size();
  for (int r=0; r<m; r++) {
    int bc = BasicVarR[r];          // col. that the basic variable is in
    if ( fabs(A[r][bc] - 1.0) > EPS) {
      double p = A[r][bc];
      for (int c=0; c<n; c++) A[r][c] /= p;
      B[r] /= p;
    }
    if (fabs(C[bc]) > EPS) {
      double p = C[bc];
      for (int c=0; c<n; c++) C[c] -= A[r][c] * p;
      obj -= B[r] * p;
    }
  }
}
bool pivoting ( vector <vector <double> > & A,
		vector <double>& B, 
		vector <double>& C,
		vector <int>& BasicVarR,   // basic variable of each row
		double & obj                // objective value
		)
{
  int m = A.size(), n = C.size();
 while (1) {
    int ev = 0;                    // id of the entering variable    
    for (ev=0; ev<n; ev++) 
      if (C[ev] < -EPS) break;

    if (ev == n) break;           // optimum reached.
          int lvr = -1;                  // leaving variable, id'ed by row
    double minRatio;
    for (int r=0; r<m; r++) {
      if (A[r][ev] > EPS) {
	if ( lvr < 0 || B[r]/A[r][ev] < minRatio ) {           
	  lvr = r; minRatio = B[r] / A[r][ev]; 
	}
      }
    }
    if (lvr < 0) return true;      // unbounded
    int lv = BasicVarR[lvr];       // leaving variable
    BasicVarR[lvr] = ev;
    double p = A[lvr][ev];
    for (int c=0; c<n; c++) A[lvr][c] /= p; B[lvr] /= p;
    for (int r=0; r<m; r++) {
      if ( r != lvr && fabs (A[r][ev]) > EPS ) {
	double p2 = A[r][ev];
	for (int c=0; c<n; c++) A[r][c] -= p2 * A[lvr][c];
	B[r] -= p2 * B[lvr];
      }
    }
    if ( fabs (C[ev]) > EPS ) {
      double p2 = C[ev];
      for (int c=0; c<n; c++) C[c] -= p2 * A[lvr][c];
      obj -= p2 * B[lvr];
    }
 if (DEBUG) {
      for (int c=0; c<n; c++) cout << C[c] << "\t"; cout << obj << endl;
      for (int r=0; r<m; r++) {
	for (int c=0; c<n; c++) cout << A[r][c] << "\t";
	cout << B[r] << endl;
      }
      cout << endl;
    }
  }
  return false;
}
void LU_solver ( vector <vector <double> > & A, // matrix A
		 vector <double>& B,           // b
		 vector <double>& X            // x
		 )
{
  int n = A.size();
  if (X.size() != n)
    X.resize (n);
vector <vector <double> > L ( n, vector<double> (n) );
  vector <vector <double> > U ( n, vector<double> (n) );
  for (int i=0; i<n; i++) 
    L[i][i] = 1.0; // diagonals of L are 1's
  copy ( A[0].begin(), A[0].end(), U[0].begin() );
   for (int k=0; k<n-1; k++) {
    for (int i=k+1; i<n; i++) { // compute the k'th column of L
      double t = A[i][k];
      for (int j=0; j<k; j++)
	t -= ( L[i][j] * U[j][k] );
      L[i][k] = t / U[k][k];
    }
    for (int j=k+1; j<n; j++) { // compute the (k+1)'s row of U
      double t = A[k+1][j];
      for (int i=0; i<k+1; i++) 
	t -= ( L[k+1][i] * U[i][j] );
      U[k+1][j] = t;
    }
  }
  for (int k=0; k<n; k++) {
    X[k] = B[k];
    for (int j=0; j<k; j++) 
      X[k] -= ( X[j] * L[k][j] );
  }
  for (int k=n-1; k>=0; k--) {
    for (int j=k+1; j<n; j++)
      X[k] -= ( X[j] * U[k][j] );
    X[k] /= U[k][k];
  }
}
int preprocess ( vector <vector <double> > & A,     // constraint matrix
		 vector <double>& B,               // right hand side
		 vector <double>& X                // unknowns
		 )
{
  int m = A.size ();                 // # of constraints
  int n = A[0].size ();              // # of variables
  vector <bool> IsRedundant (m, false);  // flags for redundant constraint
  for (int r=0; r<m; r++) {
    bool allZero = true;
    for (int c=0; c<n; c++) 
      if (fabs(A[r][c]) > EPS) { allZero = false; break; }
    if (allZero) {
      if (fabs(B[r]) > EPS) return -1;
      else IsRedundant[r] = true;
    }
  }
  for (int i=0; i<m; i++) if (!IsRedundant[i]) {
    for (int j=i+1; j<m; j++) if (!IsRedundant[j]) {
      int c;
      double ratio = 0.0;
	        for (c=0; c<n; c++) {
	if ( fabs(A[i][c]) < EPS && fabs(A[j][c]) < EPS )       // both are 0 
	  continue;
	else if ( fabs(A[i][c]) < EPS && fabs(A[j][c]) > EPS || // one is 0
	     fabs(A[i][c]) > EPS && fabs(A[j][c]) < EPS )
	  break;
	else {                                            // both are nonzero
	  if ( fabs(ratio) < EPS )
	    ratio = A[i][c] / A[j][c];
	  else {
	    if ( fabs (A[i][c]/A[j][c] - ratio) > EPS ) 
	      break;
	  }
	}
      }  if (c == n) {
	if ( fabs(B[i]) < EPS && fabs(B[j]) < EPS || 
	     fabs(B[j]) > EPS && fabs (B[i]/B[j] - ratio) < EPS )
	  IsRedundant[j] = true;
	else return -1;              // inconsistency detected
      }
    }
  }
  int r;
  for (int c=0; c<n; c++) 
     r = identity_col (A, c);
if(r==-1)
count==0;
else
count==1;

  int numRedundancies = count;
  if (numRedundancies > 0) {
    int ir = 0;                      // 1 position to the right of the new A
    for (int i=0; i<m; i++) { 
      if (!IsRedundant[i]) {
	if (ir < i) {                // overiding
	  copy (A[i].begin(), A[i].end(), A[ir].begin());
	  B[ir] = B[i];
	}
	ir++;
      }
    }
    for (int i=0; i<numRedundancies; i++) { 
      A.erase (A.end()-1);
      B.erase (B.end()-1);
    }
  }
 m -= numRedundancies;
  if (m >= n) {                      // determined or overdetermined system
    vector <vector <double> > A0 (n, vector<double> (n));
    vector <double> B0 (n);
    for (int r=0; r<n; r++) {
      copy (A[r].begin(), A[r].end(), A0[r].begin());
      B0[r] = B[r];
    }

    LU_solver (A0, B0, X);
    bool nonNegative = true;
    for (int c=0; c<n; c++)
      if (X[c] < 0) { nonNegative = false; break; }
    if (!nonNegative) 
      return -1;
    bool consistent = true;
    for (int r=n; r<m; r++) {
      double lhs = 0.0;
      for (int c=0; c<n; c++) 
	lhs += A[r][c] * X[c];
      if ( fabs (lhs - B[r]) > EPS ) { // constraint c not satisfied
	consistent = false;
	break;
      }
    }
    return (consistent ? -2 : -1);
  } 
        return numRedundancies;
}
int simplex ( const vector <vector <double> > & A,  // constraint matrix 
	      const vector <double>& B,            // right hand side
	      const vector <double>& C,            // objective vector
	      vector <double>& X,                  // unknowns
	      double & obj                          // objective value
	      )
{
  int m = A.size();                  // # of inequalities
  int n = A[0].size();               // # of variables
  if (!m || m != B.size() || n != C.size()) { 
    cout << "Wrong inputs!\n"; exit(1); 
  }
  if (X.size() != n) X.resize(n);
  fill (X.begin(), X.end(), 0);
  vector <vector <double> > A0 ( m, vector<double>(n) );
  vector <double> B0 (m);
  for (int r=0; r<m; r++) 
    copy (A[r].begin(), A[r].end(), A0[r].begin() );
  copy ( B.begin(), B.end(), B0.begin() );
  int ret_val = preprocess (A0, B0, X);
  int numRedundancies;
  if (ret_val == -1)                 // inconsistent system
    return -1;
  else if (ret_val == -2)            // solved
    return 1;
  else                               // need to run Simplex
    numRedundancies = ret_val;   

  m = A0.size ();                    // size changes after redundancy removal
 vector <bool> IsBasic (n, false);  // bit flag for basic variables
  vector <int> BasicVarR (m, -1);    // basic variable of each row
  int numBasicVar = 0;
  for (int c=0; c<n; c++) {
    int r = identity_col (A, c);

    if (r >= 0 && BasicVarR[r] < 0) {
      IsBasic[c] = true;
      BasicVarR[r] = c;
      numBasicVar++;
    }
  }
  vector <vector <double> > A2 ( m, vector<double>(n) );
  vector <double> B2 (m);
  vector <double> C2 (n);  
  for (int r=0; r<m; r++) 
    copy ( A0[r].begin(), A0[r].end(), A2[r].begin() );
  copy ( B0.begin(), B0.end(), B2.begin() );
  for (int c=0; c<n; c++) 
    C2[c] = -C[c];                   // obj. vector should be negated
  obj = 0;
  if (numBasicVar < m) {
    int n1 = n;                      // Phase I need extra dummy variables
    vector <vector <double> > A1 (m, vector<double>(n) ); 
    vector <double> B1 (m);          
    vector <double> C1 (n, 0);       // new objective vector for phase I
 for (int r=0; r<m; r++)
      copy ( A0[r].begin(), A0[r].end(), A1[r].begin() );
    copy ( B0.begin(), B0.end(), B1.begin() );  // r.h.s. is the same
    for (int i=0; i<m; i++) {
      if (BasicVarR[i] < 0) {
	for (int r=0; r<m; r++) {
	  if (r == i) A1[r].push_back (1);
	  else A1[r].push_back (0);
	}
	C1.push_back (1);            
	BasicVarR[i] = n1;
	n1++;
      }
    }
     C1.resize (n1, 1);                // Adjust sizes of objective vector
canonicalize (A1, B1, C1, BasicVarR, obj);  // convert to canonical form  
    bool unbounded = pivoting (A1, B1, C1, BasicVarR, obj);  // pivoting
    if (unbounded) {
      cout << "Unbounded Phase I!" << endl;
      exit (1);
    }    bool feasible = (fabs(obj) < EPS) ? true : false;
    if (!feasible) return 0;
    for (int r=0; r<m; r++) {
      for (int c=0; c<n; c++)
	A2[r][c] = A1[r][c];
      B2[r] = B1[r];
    }
  }
  canonicalize (A2, B2, C2, BasicVarR, obj);
  bool unbounded = pivoting (A2, B2, C2, BasicVarR, obj);
 for (int r=0; r<m; r++)              // r.h.s. is the basic solution 
    X[BasicVarR[r]] = B2[r];
  return ( unbounded ? -1 : 1 );
}
main()
{
  int m, n; 
  while (1) {
        cout << "How many constraints ? "; 
    cin >> m;
cout << "How many variables ? ";
    cin >> n;
    vector <vector <double> > A (m, vector<double>(n));
    vector <double> B (m), C(n);
    cout << "Enter the  coefficients of the " << n <<" variables in the left hand side of equality constraints:\n";
    for (int i=0; i<m; i++) for (int j=0; j<n; j++) cin >> A[i][j];
    cout << "Enter the constants on the right side of equality constraints:\n"; 
    for (int i=0; i<m; i++) cin >> B[i];
    cout << "Enter the coefficients of the" << n
	<< " variables of the objective function:\n";
    for (int i=0; i<n; i++) cin >> C[i];
    vector <double> X;
    double obj;
    cout << simplex (A, B, C, X, obj) << endl;
    cout << "\nOptimal objective value = " << obj << endl;
    cout << "\nOptimal solution: ";
    for (int i=0; i<n; i++) 
      cout << X[i] << "\t";
    cout << endl;
  }
}


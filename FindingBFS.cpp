#include<iostream>
#include<vector>
#include <cmath>
using namespace std;
vector<double> gauss(vector< vector<double> > A) 
{
    int n = A.size();

    for (int i=0; i<n; i++) {
        // Search for maximum in this column
        double maxEl = abs(A[i][i]);
        int maxRow = i;
        for (int k=i+1; k<n; k++) {
            if (abs(A[k][i]) > maxEl) {
                maxEl = abs(A[k][i]);
                maxRow = k;
            }
        }

        // Swap maximum row with current row (column by column)
        for (int k=i; k<n+1;k++) {
            double tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

        // Make all rows below this one 0 in current column
        for (int k=i+1; k<n; k++) {
            double c = -A[k][i]/A[i][i];
            for (int j=i; j<n+1; j++) {
                if (i==j) {
                    A[k][j] = 0;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
        }
    }

    // Solve equation Ax=b for an upper triangular matrix A
    vector<double> x(n);
    for (int i=n-1; i>=0; i--) {
        x[i] = A[i][n]/A[i][i];
        for (int k=i-1;k>=0; k--) {
            A[k][n] -= A[k][i] * x[i];
        }
    }
    return x;
}

void print(vector< vector<double> > A) {
    int n = A.size();
    for (int i=0; i<n; i++) {
        for (int j=0; j<n+1; j++) {
            cout << A[i][j] << "\t";
            if (j == n-1) {
                cout << "| ";
            } 
        }
        cout << "\n";
    }
    cout << endl;
}

void readInput(int r)
{
	vector<double> line(r+1,0);    
	vector< vector<double> > A(r,line);
	// Read input data
	cout<<"Now enter the coefficients of the remaining variables which are not zero\n";
	for (int i=0; i<r; i++) 
	{
		for (int j=0; j<r; j++) 
		{
			cout<<"A["<<i+1<<"]["<<j+1<<"] :";			
			cin >> A[i][j];        
		}
	}
	cout<<"Now enter the coefficients of the right hand side constants\n";
	for (int i=0; i<r; i++) 
	{
		cout<<"A["<<i+1<<"]["<<r+1<<"] :";       
		cin >> A[i][r];
	}
	//Print the Augmented Matrix
	cout<<"The augmented matrix is\n";
	print(A);
	vector<double> x(r);
	x = gauss(A);
	cout<<"The values of the non-zero variables are\n";
	//cout << "Result:\t";
	for (int i=0; i<r; i++) 
	{
		cout << x[i] << " ";
	}
	cout<<endl;
	for (int i=0; i<r; i++) 
	{
		if(x[i]<0)
		{
			cout<<"This is not the Basic Feasible Solution\n";
			break;		
		}
	}
	
	
	cout << endl;
	cout<<"==================================================\n\n";

}



void combinationUtil(int arr[], int data[], int start, int end, int index, int r);
 
// The main function that prints all combinations of size r
// in arr[] of size n. This function mainly uses combinationUtil()
void printCombination(int arr[], int n, int r)
{
    // A temporary array to store all combination one by one
    int data[r];
 
    // Print all combination using temprary array 'data[]'
    combinationUtil(arr, data, 0, n-1, 0, r);
}
 
/* arr[]  ---> Input Array
   data[] ---> Temporary array to store current combination
   start & end ---> Staring and Ending indexes in arr[]
   index  ---> Current index in data[]
   r ---> Size of a combination to be printed */
void combinationUtil(int arr[], int data[], int start, int end, int index, int r)
{
    // Current combination is ready to be printed, print it
    if (index == r)
    {
	cout<<"The coefficients of the following terms are zero\n";        
	for (int j=0; j<r; j++)
            printf("%d ", data[j]);
        printf("\n");
	// Now reading the coefficients of the remaining terms
	readInput(r);
	
        return;
    }
 
    // replace index with all possible elements. The condition
    // "end-i+1 >= r-index" makes sure that including one element
    // at index will make a combination with remaining elements
    // at remaining positions
    for (int i=start; i<=end && end-i+1 >= r-index; i++)
    {
        data[index] = arr[i];
        combinationUtil(arr, data, i+1, end, index+1, r);
    }
}
 
// Driver program to test above functions
int main()
{
    
	int r,n;	
	cout<<"Input the number of Unknowns: ";
	cin>>n;
	cout<<"Input the number of Equations: ";
	cin>>r;
	int arr[n];
	for(int i=0;i<n;i++)
		arr[i]=i+1;
	printCombination(arr, n, r);
}


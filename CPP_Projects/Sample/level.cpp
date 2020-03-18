#include <iostream>
using namespace std;

void address(int *arr, int **mxptr, int **mnptr){
int i, *max, *min;

max = &arr[0];
min = &arr[0];

for (i=0;i<5;i++){
    if (*max < arr[i])
    max = &arr[i];
    if (*min > arr[i])
    min = &arr[i];
}
*mxptr = max;
*mnptr = min;
}

int main(void){
    int *p1, *p2;
    int arr[5];
    
    for(int i =0; i < 5; i++){
    cout << "Enter Integer : ";
    cin >> arr[i];
    }
address(arr, &p1, &p2);
cout << *p1 << p1 <<endl;
cout << *p2 << p2 <<endl;
}
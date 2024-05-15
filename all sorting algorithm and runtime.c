#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
#include <math.h>

// Global variables to count the number of iterations and runtime for each sorting algorithm
int b_c = 0, i_c = 0, s_c = 0, shell_c = 0, merge_c = 0, quick_c = 0, heap_c = 0, radix_c = 0;
int counting_c = 0, bucket_c = 0, comb_c = 0, gnome_c = 0;
long long b_r = 0, i_r = 0, s_r = 0, shell_r = 0, merge_r = 0, quick_r = 0, heap_r = 0, radix_r = 0;
long long counting_r = 0, bucket_r = 0, comb_r = 0, gnome_r = 0;
// Function to print the array along with iteration count
void print_A(int array[], int n, int c) {
    printf("__________\n");
    printf("|Itr:%d|\n", c);
    for (int i = 0; i < n; i++) {
        printf("%d ", array[i]);
    }
    printf("\n");
}
// Function to print the array along with a descriptive label
void print_Array(int array[], int n, const char* data) {
    printf("%s\n", data);
    for (int i = 0; i < n; i++) {
        printf("%d ", array[i]);
    }
    printf("\n");
}
// Function to reverse the array
void rev_A(int array[], int n) {
    for (int i = 0; i < n / 2; i++) {
        int temp = array[i];
        array[i] = array[n - 1 - i];
        array[n - 1 - i] = temp;
    }
}
// Function to perform Shell Sort
void fun_shell(int array[], int n) {
    int i, j, k, tmp, c = 0;
    clock_t start = clock();
    for (i = n / 2; i > 0; i = i / 2) {
        for (j = i; j < n; j++) {
            for (k = j - i; k >= 0; k = k - i) {
                if (array[k + i] >= array[k]) {
                    break;
                } else {
                    tmp = array[k];
                    array[k] = array[k + i];
                    array[k + i] = tmp;
                    c++;
                    shell_c++;
                    print_A(array, n, c);
                }
            }
        }
    }
    clock_t end = clock();
    shell_r = (long long)(end - start);
}
// Function to perform Insertion Sort
void fun_insertion(int array[], int n) {
    int t1, j, c = 0;
    clock_t start = clock();
    for (int i = 1; i < n; i++) {
        t1 = array[i];
        j = i - 1;
        while (j >= 0 && t1 < array[j]) {
            array[j + 1] = array[j];
            --j;
        }
        array[j + 1] = t1;
        c++;
        i_c++;
        print_A(array, n, c);
    }
    clock_t end = clock();
    i_r = (long long)(end - start);
}
// Function to perform Selection Sort
void fun_selection(int array[], int n) {
    int c = 0;
    clock_t start = clock();
    for (int steps = 0; steps < n; ++steps) {
        for (int i = steps + 1; i < n; ++i) {
            if (array[steps] > array[i]) {
                c++;
                s_c++;
                int temp = array[steps];
                array[steps] = array[i];
                array[i] = temp;
                print_A(array, n, c);
            }
        }
    }
    clock_t end = clock();
    s_r = (long long)(end - start);
}
// Function to perform Bubble Sort
void fun_bubble(int array[], int n) {
    int c = 0;
    clock_t start = clock();
    for (int i = 0; i < n - 1; ++i) {
        for (int j = 0; j < n - i - 1; ++j) {
            if (array[j] > array[j + 1]) {
                c++;
                b_c++;
                int temp = array[j];
                array[j] = array[j + 1];
                array[j + 1] = temp;
                print_A(array, n, c);
            }
        }
    }
    clock_t end = clock();
    b_r = (long long)(end - start);
}
// Function to merge two subarrays in Merge Sort
void merge(int array[], int left, int mid, int right) {
    int n1 = mid - left + 1;
    int n2 = right - mid;
    int L[n1], R[n2];
    for (int i = 0; i < n1; i++) L[i] = array[left + i];
    for (int j = 0; j < n2; j++) R[j] = array[mid + 1 + j];
    int i = 0, j = 0, k = left;
    while (i < n1 && j < n2) {
        if (L[i] <= R[j]) {
            array[k] = L[i];
            i++;
        } else {
            array[k] = R[j];
            j++;
        }
        k++;
    }
    while (i < n1) {
        array[k] = L[i];
        i++;
        k++;
    }
    while (j < n2) {
        array[k] = R[j];
        j++;
        k++;
    }
}
// Function to perform Merge Sort
void mergeSort(int array[], int left, int right) {
    if (left < right) {
        int mid = left + (right - left) / 2;
        mergeSort(array, left, mid);
        mergeSort(array, mid + 1, right);
        merge(array, left, mid, right);
        merge_c++;
        print_A(array, right - left + 1, merge_c);
    }
}
// Function to perform Merge Sort (wrapper function)
void fun_merge(int array[], int n) {
    clock_t start = clock();
    mergeSort(array, 0, n - 1);
    clock_t end = clock();
    merge_r = (long long)(end - start);
}

int partition(int array[], int low, int high) {
    int pivot = array[high];
    int i = (low - 1);
    for (int j = low; j <= high - 1; j++) {
        if (array[j] < pivot) {
            i++;
            int temp = array[i];
            array[i] = array[j];
            array[j] = temp;
            quick_c++;
            print_A(array, high - low + 1, quick_c);
        }
    }
    int temp = array[i + 1];
    array[i + 1] = array[high];
    array[high] = temp;
    return (i + 1);
}
// Function to perform Quick Sort
void quickSort(int array[], int low, int high) {
    if (low < high) {
        int pi = partition(array, low, high);
        quickSort(array, low, pi - 1);
        quickSort(array, pi + 1, high);
    }
}
// Function to perform Quick Sort (wrapper function)
void fun_quick(int array[], int n) {
    clock_t start = clock();
    quickSort(array, 0, n - 1);
    clock_t end = clock();
    quick_r = (long long)(end - start);
}
// Function to heapify a subtree rooted with node i which is an index in array[]
void heapify(int array[], int n, int i) {
    int largest = i;
    int left = 2 * i + 1;
    int right = 2 * i + 2;
    if (left < n && array[left] > array[largest]) largest = left;
    if (right < n && array[right] > array[largest]) largest = right;
    if (largest != i) {
        int temp = array[i];
        array[i] = array[largest];
        array[largest] = temp;
        heap_c++;
        print_A(array, n, heap_c);
        heapify(array, n, largest);
    }
}
// Function to perform Heap Sort
void fun_heap(int array[], int n) {
    clock_t start = clock();
    for (int i = n / 2 - 1; i >= 0; i--) heapify(array, n, i);
    for (int i = n - 1; i > 0; i--) {
        int temp = array[0];
        array[0] = array[i];
        array[i] = temp;
        heapify(array, i, 0);
    }
    clock_t end = clock();
    heap_r = (long long)(end - start);
}
// Function to get the maximum element from the array
int getMax(int array[], int n) {
    int max = array[0];
    for (int i = 1; i < n; i++) if (array[i] > max) max = array[i];
    return max;
}
// Function to perform counting sort based on digits place
void countSort(int array[], int n, int exp) {
    int output[n];
    int count[10] = {0};
    for (int i = 0; i < n; i++) count[(array[i] / exp) % 10]++;
    for (int i = 1; i < 10; i++) count[i] += count[i - 1];
    for (int i = n - 1; i >= 0; i--) {
        output[count[(array[i] / exp) % 10] - 1] = array[i];
        count[(array[i] / exp) % 10]--;
    }
    for (int i = 0; i < n; i++) array[i] = output[i];
}
// Function to perform Radix sort
void fun_radix(int array[], int n) {
    int max = getMax(array, n);
    clock_t start = clock();
    for (int exp = 1; max / exp > 0; exp *= 10) {
        countSort(array, n, exp);
        radix_c++;
        print_A(array, n, radix_c);
    }
    clock_t end = clock();
    radix_r = (long long)(end - start);
}
// Function to perform Counting sort
void fun_counting(int array[], int n) {
    int max = getMax(array, n);
    int output[n], count[max + 1], i;
    memset(count, 0, sizeof(count));
    clock_t start = clock();
    for (i = 0; i < n; i++) count[array[i]]++;
    for (i = 1; i <= max; i++) count[i] += count[i - 1];
    for (i = n - 1; i >= 0; i--) {
        output[count[array[i]] - 1] = array[i];
        count[array[i]]--;
        counting_c++;
        print_A(output, n, counting_c);
    }
    for (i = 0; i < n; i++) array[i] = output[i];
    clock_t end = clock();
    counting_r = (long long)(end - start);
}
// Function to perform Bucket sort
void bucketSort(int array[], int n) {
    int max = getMax(array, n);
    int bucket_count = max / 10 + 1;
    int **buckets = (int **)malloc(bucket_count * sizeof(int *));
    int *bucket_sizes = (int *)malloc(bucket_count * sizeof(int));
    for (int i = 0; i < bucket_count; i++) {
        buckets[i] = (int *)malloc(n * sizeof(int));
        bucket_sizes[i] = 0;
    }
    clock_t start = clock();
    for (int i = 0; i < n; i++) {
        int idx = array[i] / 10;
        buckets[idx][bucket_sizes[idx]++] = array[i];
    }
    for (int i = 0; i < bucket_count; i++) {
        for (int j = 1; j < bucket_sizes[i]; j++) {
            int key = buckets[i][j];
            int k = j - 1;
            while (k >= 0 && buckets[i][k] > key) {
                buckets[i][k + 1] = buckets[i][k];
                k--;
            }
            buckets[i][k + 1] = key;
            bucket_c++;
            print_A(buckets[i], bucket_sizes[i], bucket_c);
        }
    }
    int index = 0;
    for (int i = 0; i < bucket_count; i++) {
        for (int j = 0; j < bucket_sizes[i]; j++) {
            array[index++] = buckets[i][j];
        }
    }
    clock_t end = clock();
    bucket_r = (long long)(end - start);
    for (int i = 0; i < bucket_count; i++) {
        free(buckets[i]);
    }
    free(buckets);
    free(bucket_sizes);
}
// Function to perform Comb Sort
void fun_comb(int array[], int n) {
    int gap = n, swapped = 1, c = 0;
    clock_t start = clock();
    while (gap > 1 || swapped) {
        if (gap > 1) gap = (gap * 10) / 13;
        swapped = 0;
        for (int i = 0; i < n - gap; i++) {
            if (array[i] > array[i + gap]) {
                int temp = array[i];
                array[i] = array[i + gap];
                array[i + gap] = temp;
                swapped = 1;
                c++;
                comb_c++;
                print_A(array, n, c);
            }
        }
    }
    clock_t end = clock();
    comb_r = (long long)(end - start);
}
// Function to perform Gnome Sort
void fun_gnome(int array[], int n) {
    int index = 0, c = 0;
    clock_t start = clock();
    while (index < n) {
        if (index == 0 || array[index] >= array[index - 1]) {
            index++;
        } else {
            int temp = array[index];
            array[index] = array[index - 1];
            array[index - 1] = temp;
            index--;
            c++;
            gnome_c++;
            print_A(array, n, c);
        }
    }
    clock_t end = clock();
    gnome_r = (long long)(end - start);
}
// Main function
int main() {
    int array[] = {10, 8, 6, 7, 9};
    int n = sizeof(array) / sizeof(array[0]);
    
    
    printf("Happy Coding! Enjoy!!\n\n");
    // Insertion sort
    printf("--------------------------\n\n");
    print_Array(array, n, "Insertion sort Input:");
    fun_insertion(array, n);
    print_Array(array, n, "\nInsertion sort output:");
    rev_A(array, n);

    // Selection sort
    printf("\n--------------------------\n");
    print_Array(array, n, "\nSelection sort Input:");
    fun_selection(array, n);
    print_Array(array, n, "\nSelection sort output:");
    rev_A(array, n);

    // Bubble sort
    printf("\n--------------------------\n");
    print_Array(array, n, "\nBubble sort Input:");
    fun_bubble(array, n);
    print_Array(array, n, "\nBubble sort output:");
    rev_A(array, n);

    // Shell sort
    printf("\n--------------------------\n");
    print_Array(array, n, "\nShell sort Input:");
    fun_shell(array, n);
    print_Array(array, n, "\nShell sort output:");
    rev_A(array, n);

    // Merge sort
    printf("\n--------------------------\n");
    print_Array(array, n, "\nMerge sort Input:");
    fun_merge(array, n);
    print_Array(array, n, "\nMerge sort output:");
    rev_A(array, n);

    // Quick sort
    printf("\n--------------------------\n");
    print_Array(array, n, "\nQuick sort Input:");
    fun_quick(array, n);
    print_Array(array, n, "\nQuick sort output:");
    rev_A(array, n);

    // Heap sort
    printf("\n--------------------------\n");
    print_Array(array, n, "\nHeap sort Input:");
    fun_heap(array, n);
    print_Array(array, n, "\nHeap sort output:");
    rev_A(array, n);

    // Radix sort
    printf("\n--------------------------\n");
    print_Array(array, n, "\nRadix sort Input:");
    fun_radix(array, n);
    print_Array(array, n, "\nRadix sort output:");
    rev_A(array, n);

    // Counting sort
    printf("\n--------------------------\n");
    print_Array(array, n, "\nCounting sort Input:");
    fun_counting(array, n);
    print_Array(array, n, "\nCounting sort output:");
    rev_A(array, n);

    // Bucket sort
    printf("\n--------------------------\n");
    print_Array(array, n, "\nBucket sort Input:");
    bucketSort(array, n);
    print_Array(array, n, "\nBucket sort output:");
    rev_A(array, n);

    // Comb sort
    printf("\n--------------------------\n");
    print_Array(array, n, "\nComb sort Input:");
    fun_comb(array, n);
    print_Array(array, n, "\nComb sort output:");
    rev_A(array, n);

    // Gnome sort
    printf("\n--------------------------\n");
    print_Array(array, n, "\nGnome sort Input:");
    fun_gnome(array, n);
    print_Array(array, n, "\nGnome sort output:");
    rev_A(array, n);

    double i_r_t = (double)i_r / CLOCKS_PER_SEC;
    double s_r_t = (double)s_r / CLOCKS_PER_SEC;
    double b_r_t = (double)b_r / CLOCKS_PER_SEC;
    double shell_r_t = (double)shell_r / CLOCKS_PER_SEC;
    double merge_r_t = (double)merge_r / CLOCKS_PER_SEC;
    double quick_r_t = (double)quick_r / CLOCKS_PER_SEC;
    double heap_r_t = (double)heap_r / CLOCKS_PER_SEC;
    double radix_r_t = (double)radix_r / CLOCKS_PER_SEC;
    double counting_r_t = (double)counting_r / CLOCKS_PER_SEC;
    double bucket_r_t = (double)bucket_r / CLOCKS_PER_SEC;
    double comb_r_t = (double)comb_r / CLOCKS_PER_SEC;
    double gnome_r_t = (double)gnome_r / CLOCKS_PER_SEC;

    printf("\n\n");
    printf("--------------------------\n");
    printf("Algorithm Runtime summary:\n");
    printf("--------------------------\n\n");
    printf("||Insertion\t:%d-Itr\t||Time:%.6f-Sec\t|\n", i_c, i_r_t);
    printf("\n||Selection\t:%d-Itr\t||Time:%.6f-Sec\t|\n", s_c, s_r_t);
    printf("\n||Bubble   \t:%d-Itr\t||Time:%.6f-Sec\t|\n", b_c, b_r_t);
    printf("\n||Shell    \t:%d-Itr\t||Time:%.6f-Sec\t|\n", shell_c, shell_r_t);
    printf("\n||Merge    \t:%d-Itr\t||Time:%.6f-Sec\t|\n", merge_c, merge_r_t);
    printf("\n||Quick    \t:%d-Itr\t||Time:%.6f-Sec\t|\n", quick_c, quick_r_t);
    printf("\n||Heap     \t:%d-Itr\t||Time:%.6f-Sec\t|\n", heap_c, heap_r_t);
    printf("\n||Radix    \t:%d-Itr\t||Time:%.6f-Sec\t|\n", radix_c, radix_r_t);
    printf("\n||Counting \t:%d-Itr\t||Time:%.6f-Sec\t|\n", counting_c, counting_r_t);
    printf("\n||Bucket   \t:%d-Itr\t||Time:%.6f-Sec\t|\n", bucket_c, bucket_r_t);
    printf("\n||Comb     \t:%d-Itr\t||Time:%.6f-Sec\t|\n", comb_c, comb_r_t);
    printf("\n||Gnome    \t:%d-Itr\t||Time:%.6f-Sec\t|\n", gnome_c, gnome_r_t);
    

    return 0;
}

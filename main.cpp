#pragma GCC optimize("O3")
#pragma GCC optimize("Ofast")
#include <sys/time.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <omp.h>
#include <assert.h>

using namespace std;

ofstream fout("rpt.txt");

int **needleman_wunsch(string &str1, string &str2, int match_penalty, int mismatch_penalty, int gap_penalty)
{
    int m = str1.size(), n = str2.size();

    // Initialize the 2D matrix
    int **dp = new int *[m + 1];
    for (int i = 0; i <= m; i++)
        dp[i] = new int[n + 1];

    double gflops = 6 * m * n;

    struct timeval t;

    if (m < n)
    {
        swap(m, n);
        swap(str1, str2);
    }

    cout << str1 << "  " << str2 << endl;

    // Fill the first row and column with gap penalties
#pragma omp parallel for
    for (int i = 0; i <= m; i++)
        dp[i][0] = i * gap_penalty;

#pragma omp parallel for
    for (int i = 0; i <= n; i++)
        dp[0][i] = i * gap_penalty;

    for (int i = 1; i <= n; i++)
    {
        int num_iterations = min(i, m + 1);

#pragma omp parallel for
        for (int j = 0; j < num_iterations; j++)
        {
            if (j == 0 || i == j)
                continue;

            char sequence_at_iA = str1[j - 1];
            char sequence_at_jB = str2[i - j - 1];
            int diagonal_value = dp[j - 1][i - j - 1] + (sequence_at_iA == sequence_at_jB ? match_penalty : mismatch_penalty);
            int gap_A = dp[j - 1][i - j] + gap_penalty;
            int gap_B = dp[j][i - j - 1] + gap_penalty;
            dp[j][i - j] = max(max(gap_A, gap_B), diagonal_value);
        }
    }

    for (int i = 1; i <= m; i++)
    {
        int num_iterations = m + 1 - i;

#pragma omp parallel for
        for (int j = 0; j < num_iterations; j++)
        {
            if (i == -j || n == j)
                continue;

            char sequence_at_iA = str1[i + j - 1];
            char sequence_at_jB = str2[n - j - 1];
            int diagonal_value = dp[i + j - 1][n - j - 1] + (sequence_at_iA == sequence_at_jB ? match_penalty : mismatch_penalty);
            int gap_A = dp[i + j - 1][n - j] + gap_penalty;
            int gap_B = dp[i + j][n - j - 1] + gap_penalty;
            dp[i + j][n - j] = max(max(gap_A, gap_B), diagonal_value);
        }
    }

    for (int i = 0; i <= m; i++)
    {
        for (int j = 0; j <= n; j++)
        {
            cout << dp[i][j] << " ";
            fout << dp[i][j] << " ";
        }
        cout << endl;
        fout << "\n";
    }
    fout << "\n";
    fout << "SCORE: "<<dp[m][n] << " ";

    return dp;
}

int main()
{

    string A = "AAAC";
    string B = "AGC";

    cout << A << "  " << B << endl;

    int match_penalty = 1, mismatch_penalty = -1, gap_penalty = -2;
    int **matriz = needleman_wunsch(A, B, match_penalty, mismatch_penalty, gap_penalty);

    fout.close();
}

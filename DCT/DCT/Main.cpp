#include <iostream>
#include <math.h>
#include "constants.h"
#include "Image.h"
#include "Dtc.h"


double PI = 3.14159265359;
double sqr2 = 1.41421356237;
const int N = 8;
void PrintMatrix8x8(double a[8][8])
{
    for (int i = 0; i < 8; i++)
    {
        for (int j = 0; j < 8; j++)
        {
            std::cout << a[i][j] << ", ";
        }
        std::cout << std::endl;
    }
}
void PrintMatrix8x8x8x8(double a[8][8][8][8])
{
    std::cout << "{" << std::endl;
    for (int i = 0; i < 8; i++)
    {
        std::cout << "{" << std::endl;
        for (int j = 0; j < 8; j++)
        {
            std::cout << "{" << std::endl;
            for (int k = 0; k < 8; k++)
            {
                std::cout << "{";
                for (int l = 0; l < 8; l++)
                {
                    if (l != 7)
                        std::cout << a[i][j][k][l] << ", ";
                    else
                        std::cout << a[i][j][k][l];
                }
                if (k != 7)
                    std::cout << "},";
                else
                    std::cout << "}";
            }
            if (j != 7)
                std::cout << "},";
            else
                std::cout << "}";
            std::cout << std::endl;
        }
        if (i != 7)
            std::cout << "},";
        else
            std::cout << "}";
        std::cout << std::endl;
    }
    std::cout << "}" << std::endl;
}
double input_matrix[8][8] = {
        {140,144,147,140,140,155,179,175},
        {144,152,140,147,140,148,167,179},
        {152,155,136,167,163,162,152,172},
        {168,145,156,160,152,155,136,160},
        {162,148,156,148,140,136,147,162},
        {147,167,140,155,155,140,136,162},
        {136,156,123,167,162,144,140,147},
        {148,155,136,155,152,147,147,136}
};



double CI(int i)
{
    if (i == 0)
        return 1.0 / sqrt(2);
    return 1;
}
double CI2(int i)
{
    if (i == 0)
        return 1.0 / sqrt(8);
    return 2.0 / sqrt(8);
}
double DCTxx(int u, int v)
{

    double sum = 0.0;
    for (int x = 0; x < 8; x++)
    {
        for (int y = 0; y < 8; y++)
        {
            double cosx = cos(((2.0 * double(x) + 1.0) * double(u) * PI) / (2.0 * double(8)));
            double cosy = cos(((2.0 * double(y) + 1.0) * double(v) * PI) / (2.0 * double(8)));
            double cosss = cosx * cosy * input_matrix[x][y];
            sum = sum + cosss;
        }
    }
    return sum;
}
double alfa(int i) {
    if (i == 0)
        return 1.0 / sqrt(8);
    return 0.5;
}
double output_matrix[8][8];
double output_matrix2[8][8];
void DCT(double(&input)[8][8], double(&output)[8][8])
{
    for (int u = 0; u < 8; u++)
    {
        for (int v = 0; v < 8; v++)
        {
            double sum = 0.0;
            for (int x = 0; x < 8; x++)
            {
                for (int y = 0; y < 8; y++)
                {
                    double cosx = cos((PI * (2.0 * double(x) + 1) * double(u)) / 16.0);
                    //std::cout << "cosx=" << cosx << std::endl;
                    double cosy = cos((PI * (2.0 * double(y) + 1) * double(v)) / 16.0);
                    //std::cout << "cosy=" << cosy << std::endl;
                    double pom = input[x][y] * cosx * cosy;
                    //std::cout << "pom=" << pom << std::endl;
                    sum = sum + pom;
                    //std::cout << "sum=" << sum << std::endl;
                }
            }
            output[u][v] = round(sum * alfa(u) * alfa(v));
        }
    }
}
void DCT2(double(&input)[8][8], double(&output)[8][8])
{
    for (int u = 0; u < 8; u++)
    {
        for (int v = 0; v < 8; v++)
        {
            double sum = 0.0;
            for (int x = 0; x < 8; x++)
            {
                for (int y = 0; y < 8; y++)
                {
                    sum = sum + input[x][y] * constants::cosines[u][v][x][y];
                }
            }
            output[u][v] = sum * constants::coefficients[u][v];
        }
    }
}
void IDCT(double(&input)[8][8], double(&output)[8][8])
{
    for (int x = 0; x < 8; x++)
    {
        for (int y = 0; y < 8; y++)
        {
            double sum = 0.0;
            for (int u = 0; u < 8; u++)
            {
                for (int v = 0; v < 8; v++)
                {
                    double cosx = cos((PI * (2.0 * double(x) + 1) * double(u)) / 16.0);
                    double cosy = cos((PI * (2.0 * double(y) + 1) * double(v)) / 16.0);
                    double pom = input[u][v] * cosx * cosy * alfa(u) * alfa(v);
                    sum = sum + pom;
                }
            }
            output[x][y] = sum;
        }
    }
}
void IDCT2(double(&input)[8][8], double(&output)[8][8])
{
    for (int x = 0; x < 8; x++)
    {
        for (int y = 0; y < 8; y++)
        {
            double sum = 0.0;
            for (int u = 0; u < 8; u++)
            {
                for (int v = 0; v < 8; v++)
                {
                    sum = sum + constants::cosines[u][v][x][y] * constants::coefficients[u][v] * input[u][v];
                }
            }
            output[x][y] = round(sum);
        }
    }
}
int main() {
    /*PrintMatrix8x8(input_matrix);
    DCT2(input_matrix, output_matrix);
    PrintMatrix8x8(output_matrix);
    std::cout << std::endl;
    IDCT2(output_matrix, output_matrix2);
    PrintMatrix8x8(output_matrix2);

    std::cout << constants::masterczulki << std::endl;*/

    Image test("C:\\Users\\ms\\Desktop\\obrazy\\test2.jpg");
    printf("%d, %d, %d, %d  \n", test.h,test.w, test.size, test.channels);
    for (int i = 0; i < test.w * test.channels * test.h; ++i)
    {
        if(i%3==1)
        printf("%03d, ", test.data[i]);
    }
    printf("\n\n\n ");
    int* greenArray = test.getGreenArray();

    for (int i = 0; i < (test.w * test.h); ++i)
    {
        printf("%d, ", greenArray[i]);
    }
    printf("\n\n\n ");

    //test.putGreenArrayIntoData(greenArray);

    Dtc dtc(greenArray, test.w, test.h);
    dtc.printDtc();
    dtc.writeMessage("abcd");
}
// qmc.cpp : �R���\�[�� �A�v���P�[�V�����̃G���g�� �|�C���g���`���܂��B
//


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fstream>
//#include <string.h> //�ǉ�

/* if you don't have drand48 uncomment the following two lines */
#define drand48 1.0/RAND_MAX*rand 
#define srand48 srand

#define max 250000			/* number of trials */
#define seed 68111			/* seed for number generator */


/* function returns energy of the path configuration */
double energy(double array[])
{
    int i;
    double sum = 0.;

    for (i = 0; i < 100; i++)
    {
        sum += pow(array[i + 1] - array[i], 2.0) + array[i] * array[i];
    }
    return (sum);
}


void gnuplot() {
    // PROGRA~1 ��Program Files(x86)
    FILE* gp = _popen("C:/PROGRA~1/gnuplot/bin/gnuplot.exe -persist", "w");
    if (gp == NULL) return;
    fputs("load 'D:/VCPP/VC++2022/computationPhysics/path_integral/pathIntegal_1D/pathIntegal_1D/gChart_2D.txt'\n", gp);
    fflush(gp);
    _pclose(gp);
}


int main()
{
    // path[j]��index j�����ԍ���t_j�ɑ����Bpath[j]������t_j�ł̊i�q�ʒuelement
    // prob[element]��element�Ԗڂ̊i�q�ʒu�ɗ��q�����݂���count���B
    // ���̃v���O�����ł͎��Ԋi�q���A�ʒu�i�q���͓����l��100�ɂƂ��Ă���
    //�ȉ��ł͕ϐ�element��path��prob�����̓Y�����Ɏg���Ă���̂ō������₷��
    double change, newE, oldE, path[101];
    int i, j, element, prop[101];

    double energy(double array[]);	   /* finds energy of path */

    FILE* output;                           /* save data in qmc.dat */
    //output = fopen("qmc.dat","w");
    fopen_s(&output, "D:/VCPP/VC++2022/computationPhysics/path_integral/pathIntegal_1D/pathIntegal_1D/qmc.txt", "w");

    srand48(seed);			   /* seed number generator */

    for (j = 0; j <= 100; j++) path[j] = 0.0;     /* initial path */
    for (j = 0; j <= 100; j++) prop[j] = 0;       /* initial probability */

    oldE = energy(path);                    /* find energy of path */

    for (i = 0; i < max; i++)
    {
        element = drand48() * 101;             /* pick one random element */
        change = (int)((drand48() - 0.5) * 20) / 10.0;   /* change -0.9..0.9 */

        path[element] += change;		   /* change path */

        newE = energy(path);                 /* find the new energy  */

        if ((newE > oldE) && (exp(-newE + oldE) <= drand48()))
        {
            path[element] -= change;	   /* reject */
        }

        for (j = 0; j <= 100; j++)              /* add up probabilities  */
        {
            element = path[j] * 10 + 50;
            prop[element]++;
        }
        oldE = newE;
    }

    for (i = 0; i <= 100; i++)
    {
        //fprintf(output, "%d\t%f\n", i-50, (double) prop[i]/max);
        fprintf(output, "%d %f\n", i - 50, (double)prop[i] / max);
    }
    printf("data stored in qmc.dat\n");
    fclose(output);

    gnuplot();
    return 0;

}



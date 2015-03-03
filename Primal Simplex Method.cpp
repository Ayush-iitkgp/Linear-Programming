#include <stdio.h>
#include <math.h>
#define CMAX 10
#define VMAX 10
int NC, NV, NOPTIMAL,P1,P2,XERR;
double TS[CMAX][VMAX];
void Data() {
double R1,R2;
char R;
int I,J;
printf("\n SIMPLEX METHOD\n\n");
printf(" MAXIMIZE (Y/N) ? "); scanf("%c", &R);
printf("\n NUMBER OF VARIABLES OF THE FUNCTION ? ");
scanf("%d", &NV);
printf("\n NUMBER OF CONSTRAINTS ? "); scanf("%d", &NC);
if (R == 'Y' || R=='y')
R1 = 1.0;
else
R1 = -1.0;
printf("\n INPUT COEFFICIENTS OF THE FUNCTION:\n");
for (J = 1; J<=NV; J++) {
printf(" #%d ? ", J); scanf("%lf", &R2);
TS[1][J+1] = R2 * R1;
}
printf(" Right hand side ? "); scanf("%lf", &R2);
TS[1][1] = R2 * R1;
for (I = 1; I<=NC; I++) {
printf("\n CONSTRAINT #%d:\n", I);
for (J = 1; J<=NV; J++) {
printf(" #%d ? ", J); scanf("%lf", &R2);
TS[I + 1][J + 1] = -R2;
}
printf(" Right hand side ? "); scanf("%lf", &TS[I+1][1]);
}
printf("\n\n RESULTS:\n\n");
for(J=1; J<=NV; J++) TS[0][J+1] = J;
for(I=NV+1; I<=NV+NC; I++) TS[I-NV+1][0] = I;
}
void Pivot();
void Formula();
void Optimize();
void Simplex() {
e10: Pivot();
Formula();
Optimize();
if (NOPTIMAL == 1) goto e10;
}
void Pivot() {
double RAP,V,XMAX;
int I,J;
XMAX = 0.0;
for(J=2; J<=NV+1; J++) {
if (TS[1][J] > 0.0 && TS[1][J] > XMAX) {
XMAX = TS[1][J];
P2 = J;
}
}
RAP = 999999.0;
for (I=2; I<=NC+1; I++) {
if (TS[I][P2] >= 0.0) goto e10;
V = fabs(TS[I][1] / TS[I][P2]);
if (V < RAP) {
RAP = V;
P1 = I;
}
e10:;}
V = TS[0][P2]; TS[0][P2] = TS[P1][0]; TS[P1][0] = V;
}
void Formula() {;
int I,J;
for (I=1; I<=NC+1; I++) {
if (I == P1) goto e70;
for (J=1; J<=NV+1; J++) {
if (J == P2) goto e60;
TS[I][J] -= TS[P1][J] * TS[I][P2] / TS[P1][P2];
e60:;}
e70:;}
TS[P1][P2] = 1.0 / TS[P1][P2];
for (J=1; J<=NV+1; J++) {
if (J == P2) goto e100;
TS[P1][J] *= fabs(TS[P1][P2]);
e100:;}
for (I=1; I<=NC+1; I++) {
if (I == P1) goto e110;
TS[I][P2] *= TS[P1][P2];
e110:;}
}
void Optimize() {
int I,J;
for (I=2; I<=NC+1; I++)
if (TS[I][1] < 0.0) XERR = 1;
NOPTIMAL = 0;
if (XERR == 1) return;
for (J=2; J<=NV+1; J++)
if (TS[1][J] > 0.0) NOPTIMAL = 1;
}
void Results() {
int I,J;
if (XERR == 0) goto e30;
printf(" NO SOLUTION.\n"); goto e100;
e30:for (I=1; I<=NV; I++)
for (J=2; J<=NC+1; J++) {
if (TS[J][0] != 1.0*I) goto e70;
printf(" VARIABLE #%d: %f\n", I, TS[J][1]);
e70: ;}
printf("\n ECONOMIC FUNCTION: %f\n", TS[1][1]);
e100:printf("\n");
}
int main() 
{
  Data();
  Simplex();
  Results();
  return 0;
}

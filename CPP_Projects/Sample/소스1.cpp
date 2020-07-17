#include <stdio.h>

int main() {
	
	int a = 5;
	int b = 3;

	int hap = a + b;
	int cha = a - b;
	int gop = a * b;
	int mok = a / b;
	int namuji = a % b;

	printf('%d + %d=%d/n', a, b, hap);
	printf('%d - %d=%d/n', a, b, cha);
	printf('%d * %d=%d/n', a, b, gop);
	printf('%d / %d=%d/n', a, b, mok);
	printf('%d %% %d=%d/n', a, b, namuji);

}
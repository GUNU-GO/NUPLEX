#include <stdio.h>

int main() {
	char a;
	scanf_s("%c", &a);
	printf("%d\n", a);
	printf("당신이 입력한 문자는 %c 입니다.", a);
}
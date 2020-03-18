#include <iostream>

class Point {
  int x, y;

 public:
  Point(int pos_x, int pos_y){
    x = pos_x;
    y = pos_y;
  }
};

class Geometry {
  // �� 100 ���� �����ϴ� �迭.
  Point* point_array[100];
  int point_num = 0;
 public:
  Geometry(Point **point_list);
  Geometry();

  void AddPoint(const Point &point){
    point_array[point_num]  = new Point(point);
    point_num += 1;
  }
  // ��� ���� ���� �Ÿ��� ����ϴ� �Լ� �Դϴ�.

  void PrintPoints(){
    for(int i = point_num; i>0 ;i--){
      std::cout << point_array[point_num] << std::endl;
  }
  void PrintDistance();
  // ��� ������ �մ� ������ ���� ������ ���� ������ִ� �Լ� �Դϴ�.
  // ���������� ������ �� ���� �մ� ������ �������� f(x,y) = ax+by+c = 0
  // �̶�� �� �� ������ �ٸ� �� �� (x1, y1) �� (x2, y2) �� f(x,y)=0 �� ��������
  // ���� �ٸ� �κп� ���� ������ f(x1, y1) * f(x2, y2) <= 0 �̸� �˴ϴ�.
  void PrintNumMeets();
}
};

int main(void){
Point p1(1,1);
Point p2(2,2);
Geometry geo1;
geo1.AddPoint(p1);
geo1.PrintPoints();
return 0;
}

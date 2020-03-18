#include <iostream>

class Date {
  int year_;
  int month_;  // 1 부터 12 까지.
  int day_;    // 1 부터 31 까지.

 public:
  void SetDate(int year, int month, int day);
  void AddDay(int inc);
  void AddMonth(int inc);
  void AddYear(int inc);
  void ShowDate();
};

int main(){
  Date mydate;
  mydate.SetDate(1994,6,26);
  mydate.ShowDate();
  mydate.AddDay(10);
  mydate.ShowDate();
  //mydate.AddMonth(6);
  //mydate.ShowDate();
  //mydate.AddYear(9);
  //mydate.ShowDate();

  return 0;
}

void Date::SetDate(int year, int month, int day){
  year_ = year;
  month_ = month;
  day_ = day;
}

void Date::AddDay(int inc){
  day_ += inc;
  if (day_>30)
  month_ += 1;
  day_ -= 30;
}

void Date::ShowDate(){
  std::cout << year_ << "." << month_ << "." << day_ << std::endl;
}
template<int kStepCount = 1000000>
db SimpsonIntegration(std::function<db(db)> f, db x1, db x2) {
  const db step = (x2 - x1) / kStepCount;
  const db half_step = step / 2;
  db res = 0;
  for (db x = x1; x <= x2; x += step) {
    res += (f(x) + 4 * f(x + half_step) + f(x + step));
  }
  return res * step / 6;
}

// example

db sq(db x) {
  return x * x;
}

SimpsonIntegration(sq, 0, 5);

// or

SimpsonIntegration([](db x) {
  return x * x;
}, 0, 5);


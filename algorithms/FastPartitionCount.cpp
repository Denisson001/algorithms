// Euler's pentagonal number theorem
// O(n * sqrt(n)) algo
// func[4] = 5
// 4 = 1 + 1 + 1 + 1
//     1 + 1 + 2
//     1 + 3
//     2 + 2
//     4

const int N = 3e5;
int func[N];

func[0] = func[1] = 1;
for (int i = 2; i < N; ++i) {
  for (int q = 1; ; ++q) {
    int w1 = (3 * q * q + q) / 2;
    int w2 = (3 * q * q - q) / 2;
    
    if (max(i - w1, i - w2) < 0) break;

    if (q % 2 == 1) {
      if (i - w1 >= 0) func[i] += func[i - w1];
      if (i - w2 >= 0) func[i] += func[i - w2];
    } else {
      if (i - w1 >= 0) func[i] -= func[i - w1];
      if (i - w2 >= 0) func[i] -= func[i - w2];
    }
  }
}

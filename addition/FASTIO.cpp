const int MAX_MEM = 1e8;
int mpos = 0;
char mem[MAX_MEM];
void * operator new ( size_t n ) {
    char *res = mem + mpos;
    mpos += n;
    return (void *)res;
}
void operator delete ( void * ) { }

1. /** Interface */
2.  
3. inline int readChar();
4. template <class T = int> inline T readInt();
5. template <class T> inline void writeInt( T x, char end = 0 );
6. inline void writeChar( int x );
7. inline void writeWord( const char *s );
8.  
9. /** Read */
10.  
11. static const int buf_size = 4096;
12.  
13. inline int getChar() {
14.     static char buf[buf_size];
15.     static int len = 0, pos = 0;
16.     if (pos == len)
17.         pos = 0, len = fread(buf, 1, buf_size, stdin);
18.     if (pos == len)
19.         return -1;
20.     return buf[pos++];
21. }
22.  
23. inline int readChar() {
24.     int c = getChar();
25.     while (c <= 32)
26.         c = getChar();
27.     return c;
28. }
29.  
30. template <class T>
31. inline T readInt() {
32.     int s = 1, c = readChar();
33.     T x = 0;
34.     if (c == '-')
35.         s = -1, c = getChar();
36.     while ('0' <= c && c <= '9')
37.         x = x * 10 + c - '0', c = getChar();
38.     return s == 1 ? x : -x;
39. }
40.  
41. /** Write */
42.  
43. static int write_pos = 0;
44. static char write_buf[buf_size];
45.  
46. inline void writeChar( int x ) {
47.     if (write_pos == buf_size)
48.         fwrite(write_buf, 1, buf_size, stdout), write_pos = 0;
49.     write_buf[write_pos++] = x;
50. }
51.  
52. template <class T>
53. inline void writeInt( T x, char end ) {
54.     if (x < 0)
55.         writeChar('-'), x = -x;
56.  
57.     char s[24];
58.     int n = 0;
59.     while (x || !n)
60.         s[n++] = '0' + x % 10, x /= 10;
61.     while (n--)
62.         writeChar(s[n]);
63.     if (end)
64.         writeChar(end);
65. }
66.  
67. inline void writeWord( const char *s ) {
68.     while (*s)
69.         writeChar(*s++);
70. }
71.  
72. struct Flusher {
73.     ~Flusher() {
74.         if (write_pos)
75.             fwrite(write_buf, 1, write_pos, stdout), write_pos = 0;
76.     }
77. } flusher;
78.  
79. /** Example */

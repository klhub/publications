function x = xi(x)

x1 = x(1);
x2 = x(2);
x3 = x(3);
x4 = x(4);

x  = [  x4 -x3  x2
        x3  x4 -x1
       -x2  x1  x4
       -x1 -x2 -x3 ];

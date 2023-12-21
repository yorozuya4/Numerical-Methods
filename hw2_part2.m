function y = CBE206_hw2_20190844_part2(m,N)
%section 1
A=zeros(N,N);
b=zeros(N,1);
b(N)=1;

for i=1:N
    for j=1:N
        if (j<i-m)||(j>i+m)
            A(i,j)=0;
        else
            A(i,j)=randn; 
        end
    end
end
fprintf("A = \n");
disp(A)

%section 2
tic
x1=zeros(N,1);

%Gaussian Elimination
for i=N:-1:2
    for j=i-1:-1:1
        ratio=A(j,i)/A(i,i);
        for k=1:i
            A(j,k)=A(j,k)-ratio*A(i,k);
        end
        b(j)=b(j)-ratio*b(i);
    end
end

%back substitution
for i=1:N
    remainder=0;
    for j=1:i-1
        remainder=remainder+A(i,j)*x1(j);
    end
    x1(i)=(b(i)-remainder)/A(i,i);
end
wtime1 = toc;

fprintf('after Gaussian Elimination A = \n');
disp(A)
fprintf("x1 vector from Gaussian elimination + back substitution \n");
disp(transpose(x1))
fprintf(1, "My original program took %f seconds to run.\n", wtime1);

%section 3
tic
x2=zeros(N,1);

%Gaussian Elimination
for i=N:-1:2
    if A(i,i)~=0
        for j=i-1:-1:1
            if A(j,i)~=0
                ratio=A(j,i)/A(i,i);
                for k=1:i
                    A(j,k)=A(j,k)-ratio*A(i,k);
                end
                b(j)=b(j)-ratio*b(i);
            end
        end
    end
end

%back substitution
for i=1:N
    remainder=0;
    for j=1:i-1
        remainder=remainder+A(i,j)*x2(j);
    end
    x2(i)=(b(i)-remainder)/A(i,i);
end
wtime2 = toc;

fprintf("x2 vector from Gaussian elimination + back substitution \n");
disp(transpose(x2))
fprintf(1, "My optimized program took %f seconds to run.\n", wtime2);

fprintf("sum1 = %f , sum2 = %f", sum(x1), sum(x2));

end
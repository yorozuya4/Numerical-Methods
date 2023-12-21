function y = CBE206_hw2_20190844_part1(N,A,b)
    x=zeros(N,1);

    for i=1:N-1
        %check for the suitable pivot at each column of index i
        candi_pivot=abs(A(i,i));
        if candi_pivot<1e-10
            fprintf("pivot %d is too small for index %d \n", candi_pivot,i);
            %initialize the inversibility of matrix to false because the initial pivot was too small
            inversibility=false;
            %check starting from next row(i+1:N) to determine suitable pivot
            for row=i+1:N
                if abs(A(row,i))>=1e-10
                    inversibility=true;
                    fprintf("found pivot to switch for row %d\n", row);
                    %switch row of A
                    A([i,row],:)=A([row,i],:);
                    %switch row of B accordingly
                    b([i,row])=b([row,i]);
                    break; %no need to continue checking, go for the next column index pivot
                end
            end
            % if there is no suitable pivot after looping through all rows
            if inversibility==false
                disp("The matrix doesn't have inverse!!")
                return; %stop the function
            end
        end

        %Gaussian elimination
        for j=i+1:N
            ratio=A(j,i)/A(i,i);
            for k=i:N
                A(j,k)=A(j,k)-ratio*A(i,k);
            end
            b(j)=b(j)-ratio*b(i);
        end
    end
        %back substitution
    for i=N:-1:1
        remainder=0;
        for j=i+1:N
            remainder=remainder+A(i,j)*x(j);
        end
        x(i)=(b(i)-remainder)/A(i,i);
    end
    
    disp('The final modified matrix A is');
    disp(A);
    disp('The final solution for x is');
    disp(x); 
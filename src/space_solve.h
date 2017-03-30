
void setup_stencil();

void compute_residual();



void step_stencil(){

    if(scheme_order==1)
    {
        stencil_index = new int[scheme_order+1];
        stencil_index=[0,-1];           //[j,j-1];
    } else if (scheme_order==2) {
        stencil_index = [1,0,-1];       //[j+1,j,j-1];
    }
        elseif (scheme_order==3)
                % for 3rd order fully upwind:
            %stencil_index = [0,-1,-2,-3];  %[j,j-1,j-2,j-3];
        % for 3rd order biased upwind:
            stencil_index = [1,0,-1,-2];    %[j+1,j,j-1,j-2];
        elseif (order==4)
                stencil_index = [2,1,0,-1,-2];  %[j+2,j+1,j,j-1,j-2];
        end

                S = stencil_index;
        n=length(S);

        A = zeros(n,n);

        A(1,:)=ones(1,n);
        B = zeros(n,1);
        B(2,1)=-1;

        for i=2:n

                for j=1:n
                A(i,j) = (- S(j))^(i-1) ;
        end
                end

                FD_coeff =(A\B);

        for i=1:n
                C(i) = FD_coeff(i,1) .* exp(1i.*K.*S(i));
        end

                Asd_FD = - sum(C);

    return;
}

void compute_residual(){




    return;
}

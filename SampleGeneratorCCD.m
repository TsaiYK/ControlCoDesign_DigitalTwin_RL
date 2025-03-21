function [X,X_fea_family] = SampleGeneratorCCD(num_samples,lb_x,ub_x,A,B_fun,Uc_vertex)
    % Generate the samples from the feasible set Xfea
    
    X = lhsdesign(num_samples,3);
    X = X.*repmat((ub_x-lb_x),num_samples,1)+repmat(lb_x,num_samples,1);
    count = 0;
    % The for loop is to make sure all the samples are within the feasible set, Xfea
    X_fea_family = cell(1,num_samples);
    for i = 1:num_samples
        B = B_fun(X(i,end));
        system = LTISystem('A', A, 'B', B);
        system.x.min = lb_x(1:2);
        system.x.max = ub_x(1:2);
        system.u.min = Uc_vertex(2);
        system.u.max = Uc_vertex(1);
        Xfea = system.invariantSet();
        X_fea_family{i} = Xfea;
        while ~(X(i,1:2)'<=Xfea) % if not, generate a new sample that is also in the feasible set
            X(i,:) = rand(1,3).*(ub_x-lb_x)+lb_x;
            count = count + 1;
        end
    end
end
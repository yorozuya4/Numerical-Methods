function CBE206_hw3_20190844()

    pop_size=500;
    mutation_rate=0.2;
    N=1000;
    crossover_prob=0.5;

    % Initialization of generation of chromosome s1 e1 and fitness function
    generation=[200*rand(pop_size,1),4*rand(pop_size,1)];
    fitness=zeros(pop_size,1);

    % Rank based probability function
    num_selected=pop_size*crossover_prob;
    rank_sum=sum(1:num_selected);
    prob_dist=zeros(1,num_selected);
    for i=1:num_selected
        prob_dist(i)=(num_selected+1-i)/rank_sum;
    end
    function idx=cdf_prob(prob_dist)
        total=0;
        rand_val=rand();
        for k=1:length(prob_dist)
            total=total+prob_dist(k);
            if total>=rand_val
                idx=k;
                return;
            end
        end
    end

    % Iteration begins
    for i=1:N
        for j=1:pop_size
            e1=generation(j,1);
            s1=generation(j,2);
            fitness(j)=error_energy(e1,s1);
        end
        % Selection of parent chromosome
        [~,sorted_idx]=sort(fitness,'ascend');
        parents=generation(sorted_idx(1:num_selected),:);
        % Crossover of parents and reproduction of children
        children=zeros(num_selected,2);
        for j=1:num_selected
            father_idx=cdf_prob(prob_dist);
            mother_idx=cdf_prob(prob_dist);
            while father_idx==mother_idx
                mother_idx=cdf_prob(prob_dist);
            end 
            beta1=rand();
            beta2=rand();
            children(j,1)=beta1*parents(father_idx,1)+(1-beta1)*parents(mother_idx,1);
            children(j,2)=beta2*parents(father_idx,2)+(1-beta2)*parents(mother_idx,2);
        end
        % Mutation of children
        rand_index=randi(num_selected,[num_selected,1]);
        for m=1:length(rand_index)
            idx=rand_index(m);
            if rand()<mutation_rate
                if rand()<0.5
                    children(idx,1)=200*rand();
                else
                    children(idx,2)=4*rand();
                end
            end
        end
        % Next generation
        generation=[parents;children];
    end
    
    % Find the best performing generation from all population after N generations
    [~,best_idx]=min(fitness);
    best_fitness=fitness(best_idx);
    best_e1=generation(best_idx,1);
    best_s1=generation(best_idx,2);
    fprintf('The best fitness function is %f and the best chromosome values are e1 = %f, s1 = %f.\n',best_fitness,best_e1,best_s1);
end
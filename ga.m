% The first event is to generate a number of locations for the nodes in the
% network. We will distribute them over a space 1000 km by 1000 km. The
% number of nodes is N
N = 10;  % number of nodes
S = 5;  % number of shortest links considered
delta = 3; % Average degree of all nodes
Npop = 8; % Number of chromosomes in the population. Needs to be divisible by four
mu = 0.15; % Mutation probability. Usually set between 0.1 and 0.2
maxloop = 1000; % number of loops for GA
newprob = 1;  % If newprob = 1, we start with a new set of locations. If we want to continue with the old set, we make it zero
cont = 0;   % If cont = 1, we will continue with the existing solution, to improve it
distweight = 1; % weight multiplying the total distance cost
connectweight = 1; % weight multiplying the total connection cost (no loose ends)
linknocost = 1; % weight multiplying the difference squared between the number of links and the preferred number
crossoverweight = 1; % weight multiplying the total crossover cost
maxhopweight = 0.3; % weight multiplying the max hop count cost
meanhopweight = 3; % weight multiplying the mean hop count cost

% Initialise arrays

close all;

rand

if newprob == 1
    for m = 1:N
        for n = 1:N
            distances(m,n) = 0;
        end
        x(m) = 0;
        y(m) = 0;
    end
    for n = 1:N
        x(n) = round(1000*rand);
        y(n) = round(1000*rand);
    end
end

if newprob == 1
    links(:,:) = 0;
    solution(:,:) = 0;
end

% Next we have to compute the Euclidean distance between each pair of nodes
% and store them in an array, called distances
for m = 1:N - 1
    for n = m + 1:N
        distances(m,n) = round(sqrt((x(m)-x(n))*(x(m)-x(n))+(y(m)-y(n))*(y(m)-y(n))));
        distances(n,m) = distances(m,n);
    end
end

% There are a total possible number of links equal to N*(N - 1)/2. But we
% will only consider the shortest S links connected to any node. However,
% it may be that some nodes will be connected to more than S links. The
% next section of code produces the links we will consider, the length of
% each link and the two nodes connected to each link. All of this is stored
% in the array "links".
k = 0;
for m = 1:N
    % Now we compare the distances of all the links in row m. First we put
    % the distances and the nodes of each link in the mth row of
    % "distances" into an array called "temp"
    for n = 1:N
        temp(n, 1) = distances(m,n);
        temp(n, 2) = m;
        temp(n,3) = n;
    end
    % We have a set of distances, including distance = 0 for m = n. We now
    % use the bubble sort to arrange them in increasing order
    for i = 1:N - 1
        for n = 1:N - 1
            if temp(n, 1) > temp(n + 1, 1)
                save = temp(n, 1);
                temp(n, 1) = temp(n + 1,1);
                temp(n + 1, 1) = save;
                save = temp(n, 2);
                temp(n,2) = temp(n + 1, 2);
                temp(n + 1, 2) = save;
                save = temp(n, 3);
                temp(n, 3) = temp(n + 1, 3);
                temp(n + 1, 3) = save;
            end
        end
    end
    
    % We now have to leave the first entry in temp, which will have length
    % zero, and take the next S entries to add them to the links we will
    % consider.
    for i = 2:S + 1
        add = 1;
        for j = 1:k
            if (temp(i, 2) == links(j, 2) & temp(i, 3) ==links(j, 3)) | (temp(i, 2) == links(j, 3) & temp(i, 3) ==links(j, 2))
                add = 0;
            end
        end
        if add == 1
            k = k + 1;
            links(k, 1) = temp(i, 1);
            links(k, 2) = temp(i, 2);
            links(k, 3) = temp(i, 3);
        end
    end
end

% The next step is to form the starting solution for the genetic algorithm (GA). If
% we have Npop chromosomes and k links to consider, then we have a solution
% array of Npop rows and k columns. Each element is a 1 or a 0. We only do
% this if we have decided to start a new problem (newprob = 1) or we have decided
% not to continue to improve the current generation. Otherwise,
% we continue to improve the existing solution.

if newprob == 1 | cont == 0
    for m = 1: Npop
        for n = 1:k
            solution(m, n) = round(rand);
        end
    end
end

% Here we enter the major loop. We will make it a for loop, and will allow
% it to run a fixed number of times

for loop = 1: maxloop
    % We now discard the lower half of the solution array and replace it with
    % the children of the top half. We do this by taking the lowest cost
    % solutions in pairs

    % We will use a bit mask to decide how the genes will be shared in the next
    % generation. Each pair of parents will have its own bitmask

    Npair = Npop/4;
    for m = 1: Npair
        for n = 1: k
            mask (n) = round (rand);
        end
        for n = 1:k
            solution (2*Npair + 2*m - 1, n) = solution (2*m - 1, n) * mask (n) + solution (2*m, n) * (1 - mask (n));
            solution (2*Npair + 2*m, n) = solution (2*m - 1, n) * (1 - mask (n)) + solution (2*m, n) * mask (n);
        end
    end

    % We now introduce some random changes in the genes. However, we keep the
    % best 25% of our present population unchanged

    for m = Npair + 1: Npop
        for n = 1: k
            if rand < mu
                solution (m,n) = 1 - solution (m,n);
            end
        end
    end
    
    % We now calculate the cost function for each chromosome. The cost
    % function consists of several components. The first component is simply
    % the length of the link, if it is used at all. We start by initialising
    % the cost array

    for m = 1: Npop
        for n = 2: 7
            cost (m,n) = 0;
        end
    end
    for m = 1: Npop
        for n = 1: k
            cost(m, 2) = cost(m, 2) + solution(m, n) * links(n,1) / 1000;
        end
    end

    % The next cost will be for having a solution with less than two connections
    % to a node

    for m = 1: Npop
        for node = 1: N
            connections(node) = 0;
        end
        for n = 1: k
            node = links(n, 2);
            connections(node) = connections(node) + solution(m, n);
            node = links(n, 3);
            connections(node) = connections(node) + solution(m, n);
        end
        for node = 1: N
            if connections (node) < 2
                cost (m, 3) = cost (m, 3) + 1;
            end
        end
    end

    % Now we add a cost which penalises any chromosome with a different number
    % of links from the preferred number

    for m = 1: Npop
        linkno = 0;
        for n = 1: k
            linkno = linkno + solution(m, n);
        end
        cost (m, 4) = cost (m, 4) + (linkno - N * delta/2) * (linkno - N * delta/2);
        if m == 1
            bestlinkno = linkno;
        end
    end
    
    % Now we add a cost penalty when we find two links which cross each
    % other
    
    for m = 1: Npop
        for n = 1: k-1
            if solution (m, n) == 1
                for j = n + 1: k
                    if solution (m, j) == 1
                        % We find the coordinates of the ends of the two
                        % existing links
                        xa1 = x (links (n, 2));
                        ya1 = y (links (n, 2));
                        xa2 = x (links (n, 3));
                        ya2 = y (links (n, 3));
                        xb1 = x (links (j, 2));
                        yb1 = y (links (j, 2));
                        xb2 = x (links (j, 3));
                        yb2 = y (links (j, 3));
                        % We find the slope of link 'a'
                        if xa2 == xa1
                            slope = 10000;
                        else
                            slope = (ya2 - ya1)/(xa2 - xa1);
                        end
                        % We calculate the intercepts made by link 'a' and two lines
                        % passing through the ends of the other link on the
                        % 'y' axis
                        b1 = ya1 - slope * xa1;
                        b2 = yb1 - slope * xb1;
                        b3 = yb2 - slope * xb2;
                        % Now we do the same for link 'b'
                        if xb2 == xb1
                            slope = 10000;
                        else
                            slope = (yb2 - yb1)/(xb2 - xb1);
                        end
                        b4 = yb1 - slope * xb1;
                        b5 = ya1 - slope * xa1;
                        b6 = ya2 - slope * xa2;
                        % We require that the intercepts of the lines
                        % passing through the ends of the other link be on
                        % either side of the intercept of the line through
                        % the link itself
                        if (b1 - b2) * (b1 - b3) < 0 & (b4 - b5) * (b4 - b6) < 0
                            cost (m, 5) = cost (m, 5) + 1;
                        end
                    end
                end
            end
        end
    end
    
    % Now we calculate the routes from each node to each other node and
    % find the maximum and average hop count for each chromosome
    
    for m = 1: Npop
        hopsum = 0;
        maxhops = 0;
        for source = 1: N
            % Initialise hop counts (route costs) to 100
            for dest = 1: N
                hops (dest) = 100;
                included (dest) = 0;
            end
            hops (source) = 0;
            included (source) = 1;
            % Find the links connected to source node, and set their hop
            % count (route cost) to 1
            for n = 1: k
                if solution (m, n) == 1 & links (n, 2) == source
                    hops (links (n, 3)) = 1;
                end
                if solution (m, n) == 1 & links (n, 3) == source
                    hops (links (n, 2)) = 1;
                end
            end
            
            % Start a for loop which will add all the other nodes to the
            % included set
            for j = 1: N - 1
                % First we find the node, not yey included, which has the
                % lowest hop count
                minhop = 10000;
                for i = 1: N
                    if included (i) == 0 & hops (i) < minhop
                        minhop = hops (i);
                        newnode = i;
                    end
                end
                % We now add the new node to the set of included nodes.
                % Then we update the hop count by modifying the hop count
                % of the nodes we can reach from the new node
                included (newnode) = 1;
                for n = 1:k
                    if solution (m, n) ==1 & links (n, 2) == newnode & hops (newnode) + 1 < hops (links (n, 3))
                        hops (links (n, 3)) = hops (newnode) + 1;
                    end
                    if solution (m, n) == 1 & links (n, 3) == newnode & hops (newnode) + 1 < hops (links (n, 2))
                        hops (links (n, 2)) = hops (newnode) + 1;
                    end
                end
            end
            
            % We have now found the hop count from the source node to all
            % other nodes in the network. Now we update the maximum hop
            % count and the total hop count.
            for node = 1: N
                if hops (node) > maxhops
                    maxhops = hops (node);
                end
                hopsum = hopsum + hops (node);
            end
            % The end of this loop is the end of the hop count calculations
            % for one chromosome
        end
        % We now add a cost which is the maximum hop count, and another
        % cost which is the average hop count
        cost (m, 6) = maxhops;
        cost (m, 7) = hopsum / (N * (N - 1));
    end
    % Now we add the costs together and put them in the 1st column of the array

    for m = 1: Npop
        cost (m, 1) = distweight * cost (m, 2) + connectweight * cost (m, 3) + linknocost * cost (m, 4) + crossoverweight * cost(m, 5) + maxhopweight * cost (m, 6) + meanhopweight* cost (m, 7);
    end

    % We now have to find the half of the population with the lowest total
    % costs. We can achieve this by doing a bubble sort based on the total cost

    for m = 1: Npop - 1
        for j = 1: Npop - 1
            if cost (j, 1) > cost (j + 1, 1)
                for n = 1: k
                    tempsol (n) = solution (j, n);
                    solution (j, n) = solution (j + 1, n);
                    solution (j + 1, n) = tempsol (n);
                end
                for n = 1: 7
                    tempcost (n) = cost(j, n);
                    cost (j, n) = cost (j + 1, n);
                    cost (j + 1, n) = tempcost (n);
                end
            end
        end
    end
    if mod(loop,20) == 0
        mincost (loop/20) = cost (1, 1);
    end
    bestsol(loop) = cost(1,1);
end

plot(bestsol)
hold off

cost

bestdelta = 2 * bestlinkno/N;
bestdelta
% Now we plot the links which have been included in the best three networks
for m = 1:3
    figure
    for n = 1: k
        if solution(m, n) == 1
            xx(1) = x(links(n,2));
            xx(2) = x(links(n,3));
            yy(1) = y(links(n,2));
            yy(2) = y(links(n,3));
            plot(xx, yy)
            hold on
        end
    end
    hold off
end
hold off

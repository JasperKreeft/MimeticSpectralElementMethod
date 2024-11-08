function h=Marker(varargin)

h = varargin{1};

N = length(h);

prop_name(1) = {'Marker'};
prop_name(2) = {'MarkerFaceColor'};
prop_name(3) = {'MarkerEdgeColor'};
prop_name(4) = {'LineWidth'};


prop_values = cell(N,4);

for i=1:N
    if strcmp(h(i).Marker,'none')
        prop_values(i,1) = {'o'};
    else
        prop_values(i,1) = {h(i).Marker};
    end
end

if nargin>1
    Color = varargin{2};
    for i=1:N
        prop_values(i,2) = {Color};
    end
else
    for i=1:N
        prop_values(i,2) = {h(i).Color};
    end
end
prop_values(:,3) = {'k'};
prop_values(:,4) = {1};

set(h,prop_name,prop_values)
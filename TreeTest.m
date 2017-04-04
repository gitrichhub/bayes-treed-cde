
%dat = readtable('data/treedat.csv','ReadVariableNames',false);


% Make up data
rng(12)
N = 1000;
x1 = normrnd(3,3,N,1);
x2 = normrnd(0,2,N,1);
x3 = randsample({'a','b','c','d'},N,true)';
x4 = poissrnd(1,N,1);
y = 10*x1 + normrnd(0,.5,N,1);
dat = table(x1,x2,x3,x4);

rng(252)
tic
results = TreeMCMC(y,dat,10,10,25,.05,10,.75);
toc

Treeplot(results.Trees{100})

% Example with no relation to covariates
rng(1362)
y2 = normrnd(0,1,N,1);
tic
results2 = TreeMCMC(y2,dat,10,10,25,.05,10);
toc


% Three covariates (Yabo's example)
rng(23432)
N = 1200;
x1 = unifrnd(0,2,N,1);
x2 = unifrnd(0,2,N,1);
x3 = unifrnd(0,2,N,1);
dat = table(x1,x2,x3);
y = zeros(N,1);
I = dat{:,2} < 1 & dat{:,1} < 1;
y(I) = 1;
I = dat{:,2} < 1 & dat{:,1} > 1;
y(I) = 5;
I = dat{:,2} > 1 & dat{:,3} < 1;
y(I) = 9;
I = dat{:,2} > 1 & dat{:,3} > 1;
y(I) = 13;
y = y + normrnd(0,1,N,1);
y3 = y;
dat3 = dat;

rng(283912)
tic
results3 = TreeMCMC(y3,dat3,10000,10000,25,.05,10,.75);
toc
Treeplot(results3.Trees{10000})

% Multiset
rng(1038)
tic
result3multi = TreeMCMC_Multiset(y3,dat3,10,10,2,25,.05,10,.75);
toc

a = Tree(y3,dat3,[],.05,10);
node0 = Nodes(0,[],1,2,{2,1},1:N,0);
node1 = Nodes(1,0,3,4,{1,1},[],1);
node2 = Nodes(2,0,5,6,{3,1},[],1);
node3 = Nodes(3,1,[],[],[],[],2);
node4 = Nodes(4,1,[],[],[],[],2);
node5 = Nodes(5,2,[],[],[],[],2);
node6 = Nodes(6,2,[],[],[],[],2);
a.NodeIds = 0:6;
a.Allnodes(1:7) = {node0,node1,node2,node3,node4,node5,node6};
a = descendentdata(a,0,dat3);
a = llike_termnodes(a,y3);
a.Lliketree % likelihood of true tree structure.


% Simulate a mixture of trees...


rng(13)
tic
a = Tree(y,dat,[],.5,1);
b = birth(a,y,dat);
c = prune(b,y,.75);
[d,priordraw,startcont,endcont,nchange,nchange2] = change(b,y,dat,.75)
e = birth(b,y,dat);
f = swap(e,y,dat);
g = swap(f,y,dat);
h = prune(g,y);
toc

b2 = b;
b2.Allnodes{1}.Rule{2} = b.Allnodes{1}.Splitvals{1}(end-1);
[d2,priordraw,startcont,endcont,nchange,nchange2] = change(b2,y,dat,.75)

rng(25)
a = Tree(y,dat,[],.5,1);
b = birth(a,y,dat);
c = birth(b,y,dat);
d = birth(c,y,dat);
e = birth(d,y,dat);
f = change(e,y,dat);
g = prune(f,y);
h = swap(g,y,dat);






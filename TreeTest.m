
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
results = TreeMCMC(y,dat,100,100,25,.05,10);
toc

Treeplot(results.Trees{100})

% Example with no relation to covariates
rng(1362)
y2 = normrnd(0,1,N,1);
tic
results = TreeMCMC(y2,dat,10,10,25,.05,10);
toc




rng(13)
tic
a = Tree(y,dat,[],.5,1);
b = birth(a,y,dat);
c = prune(b,y);
d = change(b,y,dat);
e = birth(b,y,dat);
f = swap(e,y,dat);
g = swap(f,y,dat);
h = prune(g,y);
toc

rng(25)
a = Tree(y,dat,[],.5,1);
b = birth(a,y,dat);
c = birth(b,y,dat);
d = birth(c,y,dat);
e = birth(d,y,dat);
f = change(e,y,dat);
g = prune(f,y);
h = swap(g,y,dat);






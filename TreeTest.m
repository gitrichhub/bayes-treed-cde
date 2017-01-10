
%dat = readtable('data/treedat.csv','ReadVariableNames',false);


% Make up data
rng(12)
N = 1000;
x1 = normrnd(3,3,N,1);
x2 = normrnd(0,2,N,1);
x3 = randsample({'a','b','c','d'},N,true)';
y = 10*x1 + normrnd(0,.5,N,1);
dat = table(x1,x2,x3);

rng(13)
tic
a = Tree(y,dat);
b = birth(a,y,dat);
c = prune(b,y);
d = change(b,y,dat);
e = birth(b,y,dat);
f = swap(e,y,dat);
g = prune(f,y);
toc


h = prune(g);

j = prune(h);

k = prune(j);

l = prune(k);

% Test change functions
m = change(f,dat);
n = change(m,dat);

% Test Swap function
o = swap(n,dat);
p = swap(o,dat);
swap(a,dat)
swap(b,dat)
swap(c,dat)


% A mix of things




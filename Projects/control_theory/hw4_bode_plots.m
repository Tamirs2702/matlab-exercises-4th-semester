%% Q1
v1 = [0.5 1 0];
v2 =  [1/250 (0.6)/50 1];
H1 = tf(5.*[0.1 1],conv(v1,v2));
subplot(1,1,1); bode(H1);
%% Q2 The correct function.
 H2d = tf(100.*conv([1 20],[1 5000]),conv(conv([1 1],[1 50]),[1 200]));
bode(H2d);

%% Q2 - The rest of the functions
 H2a = tf(100.*conv([1 1],[1 5000]),conv(conv([1 0], [1 5]),[1,6]));
 subplot(3,1,1); bode(H2a);
 H2b = tf(100,conv([1 1],[1,6]));
 subplot(3,1,2); bode(H2b);
 H2c = tf(100,conv([1 1],[1 20]));
 subplot(3,1,3);bode(H2c);
 

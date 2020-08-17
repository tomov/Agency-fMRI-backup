S = dbstack();
p = S(2).file;
disp(['--------------------- BEGIN FILE ', p, '-------------------']);
code = fileread(p);
disp(code);
disp(['--------------------- END FILE ', p, '-------------------']);

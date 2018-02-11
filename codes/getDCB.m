function [b_sat] = getDCB()
[~, ipath] = ver_chk();
bias_p1c1 = [];
f_p1c1 = strjoin({ipath, 'SolarEclipse', 'biases', 'p1c11708.dcb'}, filesep);
f_p1p2 = strjoin({ipath, 'SolarEclipse', 'biases', 'p1p21708.dcb'}, filesep);
p1c1 = fopen(f_p1c1);
bias_p1p2 = [];
p1p2 = fopen(f_p1p2);
while ~feof(p1c1)
    line = fgets(p1c1);
    if line(1) == 'G'
        l = sscanf(line, '%c %02d %f %f \\n');
        bias_p1c1 = [bias_p1c1; l(2:end)'];
    end
end
fclose(p1c1);
while ~feof(p1p2)
    line = fgets(p1p2);
    if line(1) == 'G'
        l = sscanf(line, '%c %02d %f %f \\n');
        bias_p1p2 = [bias_p1p2; l(2:end)'];
    end
end
fclose(p1p2);
b_sat = (bias_p1p2(:, 2) - bias_p1c1(:, 2)) * -2.852;
end
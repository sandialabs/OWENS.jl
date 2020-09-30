% module AirfoilRead

% export AirfoilData, readaerodyn

% import Interpolations: interpolate, Gridded, Linear, GriddedInterpolation
%import Dierckx: Spline1D, evaluate

% type AirfoilData
%     alpha::Array{Float64,1}
%     cl::Array{Float64,1}
%     cd::Array{Float64,1}
% end

% type AirfoilData
%     cl::GriddedInterpolation
%     cd::GriddedInterpolation
% end

function [alpha_rad, cl, cd] = readaerodyn(filename)
    %"""currently only reads one Reynolds number if multiple exist"""
    
    alpha = [];
    cl = [];
    cd = [];
    %
    %     open(filename) do f
    %
    %     % skip header
    %     for i = 1:13
    %         readline(f)
    %     end
    %
    %     % read until EOT
    %     while true
    %         line = readline(f);
    %         if contains(line, "EOT")
    %             break
    %         end
    %         parts = split(line);
    %         push(alpha, float(parts(1)))
    %         push(cl, float(parts(2)))
    %         push(cd, float(parts(3)))
    %     end
    
    
    fid = fopen(filename);
    % skip header
    for i = 1:13
        myfgetl(fid);
    end
    
    while true
        
        line = myfgetl(fid);
        
        if contains(line, "EOT")
            break
        end
        
        parts = split(line);
        alpha = [alpha; str2double(parts(1))];
        cl = [cl; str2double(parts(2))];
        cd = [cd; str2double(parts(3))];
    end
    fclose(fid);
    
    alpha_rad = alpha.*pi/180;
end

% af = AirfoilData(alpha*pi/180.0, cl, cd);

% 1D interpolations for now.  ignoring Re dependence (which is very minor)
% afcl = interpolate((alpha*pi/180.0,), cl, Gridded(Linear()));
% afcd = interpolate((alpha*pi/180.0,), cd, Gridded(Linear()));
% af = AirfoilData(afcl, afcd);

%     afcl = Spline1D(alpha*pi/180, cl, s=0.1);
%     afcd = Spline1D(alpha*pi/180, cd, s=0.001);
%
%     function [cl, cd] = af(alpha)
%
%         cl = evaluate(afcl, alpha);
%         cd = evaluate(afcd, alpha);
%
%         %return cl, cd
%     end

%return af%, alpha*pi/180, cl, cd


% end

% filename = "airfoils/naca0015-wt.dat";
% af, a0, cl0, cd0 = readaerodyn(filename);
% alpha = linspace(-pi, pi, 2000);
% cl, cd = af(alpha);
% using PyPlot
% figure()
% plot(alpha, cl)
% plot(a0, cl0, "o")
% figure()
% plot(alpha, cd)
% plot(a0, cd0, "o")
% show()

% C++ compatable workaround for fget1 from https://www.mathworks.com/matlabcentral/answers/461159-read-text-file-line-by-line-in-deployed-application
function line  = myfgetl(fid)
    line = '';
    c = fread(fid,1,'char=>char');
    while ~feof(fid) && ~strcmp(c, sprintf('\n'))
        if ~strcmp(c,sprintf('\r'))
            line = [line,c];
        end
        c = fread(fid,1,'char=>char');
    end
end

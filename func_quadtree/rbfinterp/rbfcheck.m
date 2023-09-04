function maxdiff = rbfcheck(options)

tic;

nodes     = options.('x');
    y     = options.('y');

s = rbfinterp(nodes, options);

maxdiff = max(abs(s-y));


% if max(abs(s-y)) > 1e-3
%     fprintf("Unfortunately, RBF Check doesn't pass. \n");
%     fprintf('max|y - yi| = %e \n', max(abs(s-y)) );
%     pause;
% end
% 
% if (strcmp(options.('Stats'),'on'))
%     fprintf('%d points were checked in %e sec\n', length(y), toc);    
% end;
% fprintf('\n');

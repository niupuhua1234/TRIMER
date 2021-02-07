function [flag] = get_glpk_flag(cmpi_flag)
% GET_GLPK_FLAG  Convert a GLPK flag to a CMPI flag

switch cmpi_flag
    case {2,5,214}
        % (integer) optimal
        flag = 2;
    case {3,110,111,210}
        % infeasible
        flag = 3;
    case {4}
        % infeasible or unbounded
        flag = 4;
    case {108,308}
        % iteration limit
        flag = 7;
    case {109,209}
        % time limit
        flag = 9;
    case {104}
        % integer solution limit
        flag = 10;
    case {105,106}
        % node limit reached
        flag = 8;
    case {316,317}
        % unstable
        flag = 13;
    otherwise
        flag = 0;
end

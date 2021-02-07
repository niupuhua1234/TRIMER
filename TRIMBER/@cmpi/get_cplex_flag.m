function [flag] = get_cplex_flag(flag)
% GET_CPLEX_FLAG  Convert a CPLEX flag to a CMPI flag

switch flag
  case {1,101,102,121}
     % (integer) optimal
     flag = 2;
  case {2,3,103}
     % infeasible
     flag = 3;
  case {4}
     % infeasible or unbounded
     flag = 4;
  case {10}
     % iteration limit
     flag = 7;
  case {11,107,108}
     % time limit
     flag = 9;
  case {104}
     % integer solution limit
     flag = 10;
  case {105,106}
     % node limit reached
     flag = 8;
  otherwise
     flag = 0;
end

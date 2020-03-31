
% Run the file
test_transient1() %BE SURE TO HAVE THE SETTINGS RIGHT IN THE TEST FILE

OLD = ('./input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff_MODAL_VERIFICATION.out');
NEW = ('./input_files_test/1_FourColumnSemi_2ndPass_15mTowerExt_NOcentStiff.out');

old = fread(fopen(OLD));
new = fread(fopen(NEW));



if ~isequal(old,new)
    fprintf('%s\n','!!!TESTS FAILED!!!')
else
    fprintf('%s\n','TESTS PASSED')
end

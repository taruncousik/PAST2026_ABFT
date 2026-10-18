function checkAnswer(name, value)
% CHECKANSWER  Compare your result with the reference answer.
%   checkAnswer('s_soi', s_soi)   prints [OK] or what looks off.
persistent ref
if isempty(ref)
    files = dir('lab*_ref.mat');
    R = load(files(1).name);
    ref = R.ref;
end
if ~isfield(ref, name)
    fprintf('[??] no reference stored for "%s"\n', name); return
end
expected = ref.(name);
if ~isequal(size(value), size(expected))
    fprintf('[CHECK] %s: size is %s, expected %s\n', name, mat2str(size(value)), mat2str(size(expected))); return
end
err = max(abs(value(:) - expected(:))) / max(max(abs(expected(:))), 1e-12);
if err < 1e-6
    fprintf('[OK]    %s\n', name);
else
    fprintf('[CHECK] %s: differs from the reference (relative error %.2g)\n', name, err);
end
end

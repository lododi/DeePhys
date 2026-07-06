function iso = yymmddToISO(tok)
% YYMMDDTOISO Convert a 6-digit YYMMDD token to an ISO 'YYYY-MM-DD' string.
%   Accepts either a numeric value (e.g. 241128) or a string/char (e.g.
%   "241128"). Used both for Excel date columns and folder-name date
%   tokens, so the two sources can be compared directly.

    if isnumeric(tok)
        s = sprintf('%06d', tok);
    else
        s = char(tok);
    end

    if isempty(regexp(s, '^\d{6}$', 'once'))
        error('yymmddToISO:badFormat', 'Expected 6-digit YYMMDD value, got "%s"', s);
    end

    iso = string(sprintf('20%s-%s-%s', s(1:2), s(3:4), s(5:6)));
end
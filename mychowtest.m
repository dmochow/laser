function [h,pValue,stat,cValue,SS] = mychowtest(Data,varargin)
%CHOWTEST Chow tests for structural change
%
% Syntax:
%
%   [h,pValue,stat,cValue] = chowtest(X,y,bp)
%   [h,pValue,stat,cValue] = chowtest(Tbl,bp)
%   [h,pValue,stat,cValue] = chowtest(...,param,val,...)
%
% Description:
%
%   Chow tests assess the stability of coefficients b in a multiple linear
%   regression model of the form y = X*b + e. Data are split at specified
%   break points in bp. Coefficients are estimated in initial subsamples,
%   then tested for compatibility with data in complementary subsamples.
%
% Input Arguments:
%
%   X - numObs-by-numPreds matrix of predictor data for a multiple linear                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
%       regression model.
%
%   y - numObs-by-1 vector of response data for a multiple linear
%       regression model.
%
%   Tbl - numObs-by-numPreds+1 tabular array of data for a multiple linear
%       regression model, with predictor data X in the first numPreds
%       columns and response data y in the last column. 
%
%   Observations with missing (NaN) values in the predictors or the
%   response are removed from the data.
%
%   bp - Scalar or vector of break points to be used in the tests. Each
%      breakpoint is the index of a specific observation in the data. Test
%      subsamples are formed from observations 1:bp and (bp+1):end.
%
%   If bp is a scalar, the number of tests, numTests, is determined by the
%   common dimension of input parameter values, and the same bp is used in
%   each test. Otherwise, the length of bp determines numTests, and
%   separate tests are run for each value in bp.
%
% Optional Input Parameter Name/Value Pairs:
%
%   NAME        VALUE
%
%   'Intercept' Logical scalar or vector of length numTests indicating
%               whether or not to add an intercept when fitting the model.
%               The default is true. The number of coefficients in the
%               model, numCoeffs, is numPreds plus 'Intercept'.
%
%   'Test'      String or cell vector of strings of length numTests
%               indicating the type of test. Values are 'breakpoint' or
%               'forecast'. A 'breakpoint' test assesses coefficient
%               equality constraints directly, using an F statistic. In
%               this case, subsamples must both have size greater than
%               numCoeffs. A 'forecast' test assesses forecast performance,
%               using a modified F statistic. In this case, only the
%               initial subsample must have size greater than numCoeffs.
%               The default is 'breakpoint'.
%
%   'Coeffs'    Logical vector or array indicating which elements of b to
%               test for equality. The default is true(1,numCoeffs), to
%               test all of b. Vector values must be of length numCoeffs.
%               Array values must be of size numTests-by-numCoeffs. If
%               'Intercept' contains mixed logical values, numCoeffs is
%               numPreds+1, and values in the first column of 'Coeffs' are
%               ignored for models without an intercept.
%
%   'Alpha'     Scalar or vector of length numTests of nominal significance
%               levels for the tests. Values must be greater than zero and
%               less than one. The default is 0.05.
%
%   'Display'   String indicating whether or not to display the results of
%               the tests in the command window. Values are 'summary' and
%               'off'. The default is 'off' for a single test and 'summary'
%               for multiple tests.
%
%   Scalar and single string parameter values, other than 'Display', are
%   expanded to the size of numTests. Vector values and 'Coeffs' arrays
%   must share a common dimension, equal to numTests.
%
%   If bp or any parameter value, other than 'Coeffs', is a row vector, all
%   outputs are row vectors.
%
% Output Arguments:
%
%   h - Vector of Boolean decisions for the tests, with length equal to
%       numTests. Hypotheses are independent of the value of 'Test':
%
%       H0: Coefficients in b selected by the 'Coeffs' parameter are
%           identical across subsamples.
%
%       H1: Some coefficient in b selected by the 'Coeffs' parameter
%           exhibits significant change across subsamples.
%
%       Values of h equal to true indicate rejection of H0 in favor of H1.
%       Values of h equal to false indicate a failure to reject H0.
%
%   pValue - Vector of p-values of the statistics, with length equal to
%            numTests.
%
%   stat - Vector of statistics, with length equal to numTests. These are
%          Chow's F statistics, given in [1].
%
%   cValue - Vector of critical values for the tests, determined by
%            'Alpha', with length equal to numTests.
%
% Notes:
%
%   o The Chow test assumes continuity of the innovations variance across
%     structural breaks. Heteroscedasticity distorts the size and power of
%     the test. It is recommended that this assumption be verified before
%     using test results for inference.
%
%   o The 'breakpoint' test is a standard F test from the analysis of
%     covariance. The 'forecast' test makes use of the standard theory of
%     prediction intervals. Chow's contribution in [1] is to place both
%     tests within the framework of general linear hypotheses, and develop
%     appropriate statistics for testing subsets of coefficients
%     (CHOWTEST's 'Coeffs' parameter).
%
%   o The 'forecast' test can be applied to cases where both subsamples
%     have size greater than numCoeffs, where a 'breakpoint' test is
%     typically used. In such cases, however, the 'forecast' test may have
%     significantly reduced power relative to a 'breakpoint' test.
%     Nevertheless, [4] suggests use of the 'forecast' test in the presence
%     of unknown specification errors.
%
%   o The 'forecast' test is based on the unbiased predictions, with zero
%     mean error, that result from stable coefficients. However, zero mean
%     forecast error does not, in general, guarantee coefficient stability.
%     Thus, 'forecast' tests are most effective in checking for structural
%     breaks, rather than model continuity [3].
%
%   o Diagnostic statistics on each subsample regression (coefficient
%     estimates, standard errors, etc.) can be obtained using the methods
%     of the LinearModel class, such as FITLM.
%
% Example:
%
%   % Test the stability of an explanatory model of U.S. real GNP using
%   % World War II as a breakpoint:
%
%   load Data_NelsonPlosser
%   span = (1915 <= dates) & (dates <= 1970);
%   Mdl = DataTable(span,[4,5,10,1]);
%   datesMdl = dates(span);
%   bp = find(datesMdl == 1945);
%   h = chowtest(Mdl,bp,'Display','summary');
%
% References:
%
%   [1] Chow, G. C. "Tests of Equality Between Sets of Coefficients in Two
%       Linear Regressions." Econometrica. Vol. 28, 1960, pp. 591–605.
%
%   [2] Fisher, F. M. "Tests of Equality Between Sets of Coefficients in
%       Two Linear Regressions: An Expository Note." Econometrica. Vol. 38,
%       1970, pp. 361-66.
%
%   [3] Rea, J. D. "Indeterminacy of the Chow Test When the Number of
%       Observations is Insufficient." Econometrica. Vol. 46, 1978, p. 229.
%
%   [4] Wilson, A. L. "When is the Chow Test UMP?" The American
%       Statistician. Vol.32, 1978, pp. 66-68.
%
% See also CUSUMTEST, RECREG, FITLM

% Copyright 2015 The MathWorks, Inc.

% Handle dataset array inputs:

if isa(Data,'dataset')
    
    try
    
        Data = dataset2table(Data);
    
    catch 
    
        error(message('econ:chowtest:DataNotConvertible'))
    
    end
    
end

% Parse inputs and set defaults:

NumericData = isnumeric(Data);

parseObj = inputParser;

parseObj.addRequired('Data',@DataCheck);
if NumericData
    parseObj.addRequired('y',@yCheck);
end
parseObj.addRequired('bp',@bpCheck);    

parseObj.addParameter('Intercept',true,@interceptCheck);
parseObj.addParameter('Test','breakpoint',@testCheck);
parseObj.addParameter('Coeffs',true); % Check in main
parseObj.addParameter('Alpha',0.05,@alphaCheck);
parseObj.addParameter('Display','off');

parseObj.parse(Data,varargin{:});

Data = parseObj.Results.Data;
if NumericData
    y = parseObj.Results.y;
end
bp = parseObj.Results.bp;

iFlag = parseObj.Results.Intercept;
whichTest = parseObj.Results.Test;
if iscell(whichTest) 
    for testIdx = 1:length(whichTest)
        whichTest{testIdx} = ...
            validatestring(whichTest{testIdx},{'breakpoint','forecast'});
    end
else
    whichTest = validatestring(whichTest,{'breakpoint','forecast'});
end
coeffs = parseObj.Results.Coeffs;
alpha = parseObj.Results.Alpha;
dispFlag = validatestring(parseObj.Results.Display,{'summary','off'});

% Check parameter dimensions, expand scalars and single strings:

[numTests,rowOutput,bp,iFlag,whichTest,alpha] = ...
    sizeCheck(bp,iFlag,whichTest,alpha);

% Forgive nonlogical parameter values:

if ~islogical(iFlag)
    
    try
        
        iFlag = logical(iFlag);
        
    catch
        
        error(message('econ:chowtest:InterceptNotConvertible'))
        
    end
    
end

if ~isempty(coeffs) && ~islogical(coeffs)
    
    try
        
        coeffs = logical(coeffs);
        
    catch
        
        error(message('econ:chowtest:CoeffsNotConvertible'))
        
    end
    
end

% Get X, y:

if isnumeric(Data)
    
    X = Data;
    
    if length(y) ~= size(X,1)
             
    	error(message('econ:chowtest:ResponseVectorWrongSize'))
          
    else 
    
        y = y(:);
        
    end
    
else
    
    X = table2array(Data(:,1:(end-1)));
    y = table2array(Data(:,end));
    
end

% Remove observations with missing values:

D = [X,y];
D(any(isnan(D),2),:) = [];
X = D(:,1:end-1);
y = D(:,end);
[numObs,numPreds] = size(X); % Chow [1]: numObs = n + m

% Pad number of coefficients, if intercepts:

pad = any(iFlag);
numCoeffsPadded = numPreds + pad;

% Check 'coeffs' parameter:

if any(strcmp(parseObj.UsingDefaults,'Coeffs'))
    
    coeffs = true(numTests,numCoeffsPadded);
    
elseif isvector(coeffs)
    
    if length(coeffs) ~= numCoeffsPadded

        error(message('econ:chowtest:CoeffsVectorWrongSize'))
  
    else
    
	coeffs = coeffs(:)';
	coeffs = repmat(coeffs,numTests,1);
    
    end

elseif ((numTests > 1) && (size(coeffs,1) ~= numTests))...
        || (size(coeffs,2) ~= numCoeffsPadded)
        
	error(message('econ:chowtest:CoeffsArrayWrongSize'))
    
end

% Scalar expand single tests when size(coeffs,1) > 1

if (numTests == 1) && (size(coeffs,1) > 1)
        
	numTests = size(coeffs,1);
    
    bp = bp(ones(numTests,1));
    iFlag = iFlag(ones(numTests,1));
    whichTest = whichTest(ones(numTests,1));
    alpha = alpha(ones(numTests,1));
        
end

% Set 'Display' parameter:

if any(strcmp(parseObj.UsingDefaults,'Display'))
    
    if numTests == 1
        
        dispFlag = 'off';
        
    else
        
        dispFlag = 'summary';
        
    end
    
end

needDisplay = strcmp(dispFlag,'summary');

% Preallocate outputs:

h = false(numTests,1);
pValue = zeros(numTests,1);
stat = zeros(numTests,1);

needCValue = (nargout >= 4);
if needCValue
    cValue = zeros(numTests,1);
end

% Perform the tests:

for t = 1:numTests
    
    testType = whichTest{t};
    testIntercept = iFlag(t);
    numCoeffs = numPreds+testIntercept; % Chow [1]: p
    testBP = bp(t); % Chow [1]: n
    numObsRight = numObs-testBP; % Chow [1]: m
    testAlpha = alpha(t);
    
    % Check for sufficient observations:
    
    if ~(numObs > numCoeffs)
        
         error(message('econ:chowtest:TooFewObservations'))
    
    elseif ~(testBP > numCoeffs)
        
        error(message('econ:chowtest:TooFewObservationsLeftSide'))
        
    elseif strcmp(testType,'breakpoint') && (~(numObsRight > numCoeffs))
        
       error(message('econ:chowtest:TooFewObservationsRightSideForBreakpointTest')) 
        
    end
    
	% Prepare data:
    
    if testIntercept
        
        XTest = [ones(numObs,1),X];
        Idx1 = 1;
        
    else % Ignore value in pad column of 'Coeffs'
        
        XTest = X;
        Idx1 = pad+1;
        
    end
    
    testCoeffs = coeffs(t,Idx1:end);
    Z = XTest(:,testCoeffs);
    W = XTest(:,~testCoeffs);
    q = size(Z,2); % Chow [1]: q
    r = size(W,2); % Chow [1]: p-q
    
    if q == 0
        
       error(message('econ:chowtest:NoCoeffsToTest'))
        
    end
    
	if ~(numObsRight > r) % Assure forecast test dof1 > 0
        
        error(message('econ:chowtest:TooFewObservationsRightSideForConstraints'))
        
	end
    
	% Break data: 
    
    Z1 = Z(1:testBP,:);
    Z2 = Z((testBP+1):end,:);
    
	W1 = W(1:testBP,:);
	W2 = W((testBP+1):end,:);
    
	y1 = y(1:testBP);
	y2 = y((testBP+1):end);

    % Unrestricted regression:

    XU = [Z1 zeros(testBP,q) W1 zeros(testBP,r); ...
          zeros(numObsRight,q) Z2 zeros(numObsRight,r) W2];
    
    % If r < numObsRight < numCoeffs (forecast test only), XU is rank
    % deficient. Full rank 2*numCoeffs is reduced by numCoeffs-numObsRight.
    % In these cases, MLDIVIDE (backslash) returns a nonunique, basic
    % solution for bU, setting elements of c2 and (if r > 0) d2 equal to 0.
    % These 0 elements are not used in the forecast test statistic.
    % Moreover, the sum of squares of residuals, used to construct the F
    % statistic, is unaffected by the choice of bU.

    % Supress rank-deficiency warning:
    
    warnStruct = warning('off','MATLAB:rankDeficientMatrix');
	bU = XU\y
    warning(warnStruct);
    
    c1 = bU(1:q);
    c2 = bU((q+1):(2*q));
    
    if r > 0
    
        d1 = bU((2*q+1):(2*q+r));
        d2 = bU((2*q+r+1):end);
        
    else
        
        d1 = zeros(0,1);
        d2 = zeros(0,1);
        
    end
    
    % Restricted regression:
    
    XR = [Z1 W1 zeros(testBP,r); ...
          Z2 zeros(numObsRight,r) W2];
    bR = XR\y

    c0 = bR(1:q);
    
	if r > 0
        
        d10 = bR((q+1):(q+r));
        d20 = bR((q+r+1):end);
        
	else
        
        d10 = zeros(0,1);
        d20 = zeros(0,1);
        
	end
    
    % Compute the test statistic:
    
    v1 = Z1*c1+W1*d1-Z1*c0-W1*d10;
    ss1 = v1'*v1;
    
    v3 = y1-Z1*c1-W1*d1;
    ss3 = v3'*v3;

	if strcmp(testType,'breakpoint')
        
        v2 = Z2*c2+W2*d2-Z2*c0-W2*d20;
        ss2 = v2'*v2;
        
        v4 = y2-Z2*c2-W2*d2;
        ss4 = v4'*v4;
        
        dof1 = q;
        dof2 = numObs-2*numCoeffs;
        testStat = ((ss1+ss2)/(ss3+ss4))*(dof2/dof1);
        SS{1}=ss1; SS{2}=ss2; SS{3}=ss3; SS{4}=ss4; 
        
    else
        
        v2 = y2-Z2*c0-W2*d20;
        ss2 = v2'*v2;
        
        dof1 = numObsRight-numCoeffs+q;
        dof2 = testBP-numCoeffs;
        testStat = ((ss1+ss2)/ss3)*(dof2/dof1);
        SS{1}=ss1; SS{2}=ss2; SS{3}=ss3;
	end
    
    stat(t) = testStat;

    % Compute p-value:

    testPValue = fcdf(testStat,dof1,dof2,'upper');
    pValue(t) = testPValue;
    
    % Perform the test:

    testDecision = (testAlpha >= testPValue);
    h(t) = testDecision;

    % Compute critical value, if requested:

    if needCValue || needDisplay

       testCValue = finv(1-testAlpha,dof1,dof2);
       cValue(t) = testCValue;

    end
    
    % Update display:
    
    if needDisplay

        if all(testCoeffs)

            testCoeffString = 'All';

        else

           testCoeffString = num2str(testCoeffs);

        end

        if testDecision

            testDecisionString = 'Reject coefficient stability';

        else

            testDecisionString = 'Fail to reject coefficient stability';

        end

        if t == 1

            fprintf('\nRESULTS SUMMARY\n\n')

        end

        fprintf('***************')
        fprintf('\nTest %d\n',t)

        fprintf('\nSample size: %u',numObs)
        fprintf('\nBreakpoint: %u\n',testBP)

        fprintf('\nTest type: %s',testType)
        fprintf('\nCoefficients tested: %s\n',testCoeffString)

        fprintf('\nStatistic: %.4f',testStat)
        fprintf('\nCritical value: %.4f\n',testCValue)

        fprintf('\nP value: %.4f',testPValue)
        fprintf('\nSignificance level: %.4f\n',testAlpha)

        fprintf('\nDecision: %s\n\n',testDecisionString)
    
    end

end

% Orient output:

if rowOutput

   h = h';
   pValue = pValue';
   stat = stat';

   if needCValue

      cValue = cValue';

   end

end

%-------------------------------------------------------------------------
% Check input Data
function OK = DataCheck(Data)

if ~isnumeric(Data) && ~isa(Data,'table') 

    error(message('econ:chowtest:DataFormatUnsupported'))
    
elseif isa(Data,'table') && (size(Data,2) < 2)

    error(message('econ:chowtest:DataWithoutResponse'))

else
    
    OK = true;

end

%-------------------------------------------------------------------------
% Check input y
function OK = yCheck(y)
    
if ~isnumeric(y) || ~isvector(y)

    error(message('econ:chowtest:ResponseNonNumericVector'))
      
else

    OK = true;

end

%-------------------------------------------------------------------------
% Check input bp
function OK = bpCheck(bp)

if ~isnumeric(bp) || ~isvector(bp)

    error(message('econ:chowtest:BreakpointsNonNumericVector'))

elseif any(mod(bp,1) ~= 0) || any(bp < 1)

    error(message('econ:chowtest:BreakpointsOutOfRange'))

else

    OK = true;

end

%-------------------------------------------------------------------------
% Check value of 'Intercept' parameter
function OK = interceptCheck(iFlag)

if ~(isnumeric(iFlag) || islogical(iFlag)) || ~isvector(iFlag)
    
    error(message('econ:chowtest:InterceptNonVector'))

else

    OK = true;

end

%-------------------------------------------------------------------------
% Check value of 'Test' parameter
function OK = testCheck(testType)
    
if ~isvector(testType) || ...
   isnumeric(testType) || ...
   (iscell(testType) && any(cellfun(@isnumeric,testType)))

    error(message('econ:chowtest:TypeOfTestNonStringVector'))

else

    OK = true;

end

%-------------------------------------------------------------------------
% Check value of 'Alpha' parameter
function OK = alphaCheck(alpha)
    
if ~isnumeric(alpha) || ~isvector(alpha)

    error(message('econ:chowtest:AlphaNonNumericVector'))
    
elseif any(alpha <= 0) || any(alpha >= 1)
        
	error(message('econ:chowtest:AlphaOutOfRange'))

else

    OK = true;

end

%-------------------------------------------------------------------------
% Check parameter dimensions, expand scalars and single strings:
function [numTests,rowOutput,varargout] = sizeCheck(varargin)

numCheck = length(varargin);
varargout = cell(1,numCheck);

checkLengths = zeros(numCheck,1);

% Find parameter lengths:
for i = 1:numCheck
    
    ivar = varargin{i};
    
    if isnumeric(ivar) || islogical(ivar) || iscell(ivar)
        
        checkLengths(i) = length(ivar);
        
    else
        
        checkLengths(i) = 1; % Single string
        varargin{i} = varargin(i); % Convert to cell
        
    end
    
end

% Check nonscalar parameter lengths for commensurability, set numTests:

nonScalarIdx = find(checkLengths > 1);

if ~isempty(nonScalarIdx)
    
    ns1 = nonScalarIdx(1);
    
    if any(checkLengths(nonScalarIdx) ~= checkLengths(ns1))
        
        error(message('econ:chowtest:ParameterSizeMismatch'))

    end
    
    numTests = checkLengths(ns1);
    
else
    
    numTests = 1;

end

% Set row output:

[~,n] = cellfun(@size,varargin);
rowOutput = any(n > 1);

% Assign output:

for i = 1:numCheck

    if checkLengths(i) == 1

        varargout{i} = repmat(varargin{i},numTests,1);

    else

        varargout{i} = varargin{i};

    end

end
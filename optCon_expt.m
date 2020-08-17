function EXPT = optCon_expt(local)

    % something weird going on with exptdir and not being able to access
    % local or ncf directory - needs to be fixed!
    
    
    %creates EXPT structure for CCNL fMRI processing pipeline
    %
    % USAGE: EXPT = optCon_expt()
    %
    % INPUTS:
    %   local (optional) - true if file path is on local computer, false if on NCF
    %
    % OUTPUTS:
    %   EXPT - experiment structure with fields
    %          .TR - repetition time
    %          .create_multi - function handle for creating multi structure
    %          .modeldir - where to put model results
    %          .subject(i).datadir - directory for subject data
    %          .subject(i).functional - .nii files for runs
    %          .subject(i).structural - .nii for structural scan
    %
    % Cody Kommers, July 2016
    % Momchil Tomov, Nov 2016
    % Hayley Dorfman, Aug 2018
    
    
    print_code
    
    
    % set default parameters
    %
     if nargin < 1
         disp(nargin >1)
         [~, name] = system('hostname');
         disp(system('hostname'))
      if ~isempty(strfind(name, 'Hayley'))    
         % err on the side of falsely thinking it's NCF. Because locally
         % you will catch that mistake immediatley. On NCF, you will
         % catch it after you're already sent 100 jobs and they all
         % fail 2 days later...
         %
         local = true;
       else
         local = false;
        end
     end


    if local
        disp('local')
        %exptdir = '/Users/hayley/Dropbox (Personal)/studies/OptControl/UCR_fMRI/analyses/imaging/';
        exptdir = '/Users/hayleydorfman/Dropbox/studies/OptControl/UCR_fMRI/analyses/imaging/';%new comp

    else
        exptdir = '/ncf/gershman/Lab/Hayley/'; % on CBS central server
        disp('nope')
    end
    
    
    [subjdirs, nRuns, goodRuns, goodSubjects, subj_original_indices] = optCon_getSubjectsDirsAndRuns();
    %subjdirs
    %nRuns
    
    for subj = 1:length(subjdirs)
        %subjdir = [exptdir, 'subjects_ants/', subjdirs{subj}, '/']; %for
        %ants
        subjdir = [exptdir, 'subjects/', subjdirs{subj}, '/']; 
      
        EXPT.subject(subj).datadir = [subjdir, 'preproc'];
        
        EXPT.subject(subj).structural = 'struct.nii';
        
        i = 1;
        for run = 1:length(goodRuns{subj})
            if(goodRuns{subj}(run))
                EXPT.subject(subj).functional{i} = ['run',sprintf('%03d',run),'.nii'];
                i = i + 1;
            end 
        end
        disp(EXPT.subject(subj));
    end
    
    % TR repetition time
    EXPT.TR = 2; %seconds
    % Function handle to create subject multi structure
    EXPT.create_multi = @optCon_create_multi;
    % Where you want model output data to live
    EXPT.modeldir = [exptdir, 'glmOutput'];
    EXPT.nTR = 215; %number of TRs in a run
    EXPT.run_duration = EXPT.nTR * EXPT.TR; %run duration calculation
    
    % Where the data live, but not sure which data
    EXPT.datadir = [exptdir, 'testOutput'];


end

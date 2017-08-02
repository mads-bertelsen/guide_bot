function AddProject(action, project_dir)
switch action
    case 'new'
        % No project list, so make one from scratch
        allprojects{1} = 2;
        allprojects{2} = project_dir;
        save('projectlist.mat', 'allprojects')
    case 'add'
        % Check if project list already exists
        guidebotContents = dir;
        if sum(ismember({guidebotContents.name},'projectlist.mat')) == 1
            load('projectlist.mat')
            % Make sure a concurrent instrument run has not already added the project
            if ~strcmp(allprojects{end}, project_dir)
                allprojects{end+1} = project_dir;
                save('projectlist.mat', 'allprojects')
            end
        else % No project list, so make one from scratch
            allprojects{1} = 2;
            allprojects{2} = project_dir;
            save('projectlist.mat', 'allprojects')
        end
            
    case 'single'
        % Do nothing
    otherwise
        disp(['ERROR: ' action ' is and invalid name for options.projectlist. Valid options are the string ''new'', ''single'', or ''add''']);
end
    

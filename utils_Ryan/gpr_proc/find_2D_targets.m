function Targets = find_2D_targets(S, I, target_threshold, max_X, max_Z, target_MPP, target_MPS)
    MAX_BOX_H = round(max_Z/S.pixel_width);
    MAX_BOX_W = round(max_X/S.pixel_width);
    
    % look for peaks in the target
    In = normalize(I);
    I_dB = 10.*log10(In);
    
    % Algorithm to find "peaks"....
    MPP = target_MPP;
    MPS = round(target_MPS/S.pixel_width);
    ILM2        = islocalmax(I_dB,2,'MinProminence',MPP,'MinSeparation',MPS);
    ILM1        = islocalmax(I_dB,1,'MinProminence',MPP,'MinSeparation',MPS);
    ILM         = ILM2 & ILM1;
    inds        = find(ILM);
    [targ_rows, targ_cols] = ind2sub(size(I),inds);
    Ntargs      = numel(targ_rows);
    to_delete   = zeros(Ntargs,1,'logical');
    for ii = 1:Ntargs
        pk_ii = [targ_rows(ii),targ_cols(ii),0];
        for jj = 1:Ntargs
            if(jj~=ii && ~to_delete(jj))
                if(radial_distance(pk_ii,[targ_rows(jj),targ_cols(ii),0]) < MPS)
                    if(I(sub2ind(size(I),targ_rows(ii),targ_cols(ii))) > I(sub2ind(size(I),targ_rows(jj),targ_cols(jj))))
                        to_delete(jj)=true;
                    end
                end
            end
        end
    end
    targ_rows(to_delete)=[];
    targ_cols(to_delete)=[];
    Ntargs = numel(targ_rows);
    
    T = Target2D.empty(Ntargs,0);
    for ii = 1:Ntargs
        % find the local minimum for target ii
        ind_guess_x = targ_cols(ii);
        ind_guess_z = targ_rows(ii);
        I_th = I_dB > I_dB(ind_guess_z, ind_guess_x) + target_threshold;

        % draw the smallest possible box encapsulating all connected
        % elements
        topflag = 1; botflag = 1; leftflag = 1; rightflag = 1;
        box_top   = max(1,         ind_guess_z - 1); 
        box_bot   = min(size(I,1), ind_guess_z + 1);
        box_left  = max(1,         ind_guess_x - 1);
        box_right = min(size(I,2), ind_guess_x + 1);
        
        while(topflag || botflag || leftflag || rightflag)
            topflag = 0; botflag = 0; leftflag = 0; rightflag = 0;
            
            %calculate the peeks
            box_top_peek   = max(1,         box_top   - 1);
            box_bot_peek   = min(size(I,1), box_bot   + 1);
            box_left_peek  = max(1,         box_left  - 1);
            box_right_peek = min(size(I,2), box_right + 1);
            
            %check if top or bot has to be incremented (incl. diag)
            topflag = (sum(I_th(box_top_peek,box_left_peek:box_right_peek)) > 0);
            botflag = (sum(I_th(box_bot_peek,box_left_peek:box_right_peek)) > 0);
            
            %check if left or right has to be incrememnted(incl. diag)
            leftflag  = (sum(I_th(box_top_peek:box_bot_peek,box_left_peek)) > 0);
            rightflag = (sum(I_th(box_top_peek:box_bot_peek,box_right_peek)) > 0);
            
            
            if(abs(box_right - box_left) > MAX_BOX_W)
               rightflag = 0;
               leftflag  = 0;
            end
            if(abs(box_top - box_bot) > MAX_BOX_H)
               topflag = 0;
               botflag = 0;
            end
            
            if(topflag)
                box_top = box_top_peek;
                if(box_top == 1)
                    topflag = 0;
                end
            end
            if(botflag)
                box_bot = box_bot_peek;
                if(box_bot == size(I,1))
                   botflag = 0; 
                end
            end
            if(leftflag)
                box_left = box_left_peek;
                if(box_left == 1)
                   leftflag = 0; 
                end
            end
            if(rightflag)
                box_right = box_right_peek;
                if(box_right == size(I,2))
                   rightflag = 0; 
                end
            end
        end

        % at this point...
        T(ii) = Target2D(S, I, S.xC(box_left), S.zC(box_top), abs(S.xC(box_right)-S.xC(box_left)), abs(S.zC(box_bot)-S.zC(box_top)));
    end
    
    for ii = 1:Ntargs
        % look at the average energy in the target vs. average energy in
        % the non-target.
        I_ref = I;
        [~,x_ind] = min(abs(S.xC - T(ii).upper_left_x));
        [~,z_ind] = min(abs(S.zC - T(ii).upper_left_z));
        I_ref(z_ind + (0:T(ii).Nz-1), x_ind + (0:T(ii).Nx-1)) = nan;
        I_ref = I_ref(:);
        I_ref(isnan(I_ref))=[];
        
        T(ii).SNR = 10*log10(T(ii).peak / std(I_ref));
    end
    
    Targets = T([T.SNR] > 6);
    
end
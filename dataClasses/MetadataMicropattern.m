classdef MetadataMicropattern < Metadata
    % metadata with additional properties for Andor IQ output
    
    % ---------------------
    % Idse Heemskerk, 2016
    % ---------------------
    
    properties
        
        colRadiiMicron
        colRadiiPixel
        colMargin
    end
    
    methods
        
        function this = MetadataMicropattern(filename,skipPixSize)
                if ~exist('skipPixSize','var')
                    skipPixSize = false;
                end
          
                this = this@Metadata(filename,skipPixSize);
            
            
            % default margin outside colony to process, in pixels
            this.colMargin = 10;
        end
    end
end
% Monochromator reflectivity depends not only on the material, but also on
% the mosaic and thickness. A McStas refectivity table only provides data for a
% material at a specific mosaic and material thickness. This function
% calculates a prefactor that can be used to adjust the table values to
% match a different mosaic and thickness.

% When a reflectivity table is employed, McStas calculates the reflectivity
% at a given neutron wavevector as:

% R(wavevector) = r0*R_table(wavevector)

% where R_table is the reflectivity file specified with the monochromator
% option 'reflect' and r0 is the constant reflectivity specified with the
% monochromator option 'r0'; typically just set to 1 when a table is 
% employed. This function calculates what r0 should be in order to 
% properely adjust the table values for a different mosaic and/or
% thickness.

%%%%%%%%%%% INPUTS %%%%%%%%%%
% MonoMosaic:                           Monochromator mosaic selected by user.
% MonoThickness:                        Monochromator thickness selected by user.
% TableMosaic:                          Monochromator mosaic of the reflectivity table.
% TableThickness:                       Monochromator thickness of the reflectivity table.
% AverageTableReflectivity [0 1]:       Average table reflectivity over the range of wavelengths of interest (see note below.)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%% OUTPUTS %%%%%%%%%%
% r0:                                   Prefactor to adjust a reflectivity table to match a different thickness/mosaic.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% NOTE: Technically, each value in a reflectivity table should be
% multiiplied by a different value in order to properly adjust for a new
% mosaic/thicknes. However, if the reflectivity does not have much
% variation, then a single r0 value can be used based on the average table
% reflectivity over the range of wavelengths of interest. This function is
% written so that the AverageTableReflectivity can be a vector containing
% all the reflectivity values in a table, it will then calculate the proper
% prefactor for each table value. It is up to the user to generate the new
% table though. You would then use this new table and revert back to
% setting r0 = 1.

% NOTE: MonoMosaic, TableMosaic and MonoThickness, TableThickness can use
% any units, as long as they are consistent.

% NOTE: This calculation follows from the peak reflectivity equation reported in
% chapter 3.2 of Shirane, Shapiro, and Tranquada


% EXAMPLE:
% A user has a PG002 monochromator with a thickness of 1.5mm and a mosaic of 45min.
% The HOPG.rfl table included with McStas was taken using a monochromator
% with 2mm thickness and 30min mosaic. Neutrons with wavevector (energy)
% over the range of 1.098AA-1 (2.5meV) to 2.691AA-1 (15meV) will be studied.
% The average reflectivity over this range is ~0.78 (based on the average of the two
% endpoints in the table.) Thus, the renormalization constant r0 is:

% r0 = AdjustReflectivityTable(45, 1.5, 30, 2, 0.78)

% I am actually guessing the mosaic/thickness in the HOPG.rfl table, I have
% contacted the McStast team to get the true values, if you are reading this
% then that means I have not recieved an answer. In which case, if you know the
% mosaic/thickness used then please let me know!

% Leland Harriger, ORNL
% July 2017

function r0 = AdjustReflectivityTable(MonoMosaic, MonoThickness, TableMosaic, TableThickness, AverageTableReflectivity)
% Redefine inputs
m = MonoMosaic;
t = MonoThickness;
m0 = TableMosaic;
t0 = TableThickness;
R0 = AverageTableReflectivity;

% Do Calculation
r0 = m0*t/((1-R0)*m*t0 + R0*m0*t);









// (c) Copyright BCL @ Vanderbilt University 2014
// (c) BCL Homepage: http://www.meilerlab.org/bclcommons
// (c) BCL Code Repository: https://github.com/BCLCommons/bcl
// (c)
// (c) The BioChemical Library (BCL) was originally developed by contributing members of the Meiler Lab @ Vanderbilt University.
// (c)
// (c) The BCL is now made available as an open-source software package distributed under the permissive MIT license,
// (c) developed and maintained by the Meiler Lab at Vanderbilt University and contributing members of the BCL Commons.
// (c)
// (c) External code contributions to the BCL are welcome. Please visit the BCL Commons GitHub page for information on how you can contribute.
// (c)
// (c) This file is part of the BCL software suite and is made available under the MIT license.
// (c)

// initialize the static initialization fiasco finder, if macro ENABLE_FIASCO_FINDER is defined
#include "util/bcl_util_static_initialization_fiasco_finder.h"
BCL_StaticInitializationFiascoFinder

// include header of this class
#include "app/bcl_app_apps.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_aa_neighbor_count.h"
#include "assemble/bcl_assemble_aa_neighbor_list_container.h"
#include "assemble/bcl_assemble_aa_neighbor_vector.h"
#include "assemble/bcl_assemble_protein_ensemble.h"
#include "command/bcl_command_app_default_flags.h"
#include "command/bcl_command_command.h"
#include "command/bcl_command_flag_dynamic.h"
#include "command/bcl_command_flag_static.h"
#include "command/bcl_command_parameter_check_allowed.h"
#include "command/bcl_command_parameter_check_enumerate.h"
#include "command/bcl_command_parameter_check_ranged.h"
#include "math/bcl_math_gnuplot_heatmap.h"
#include "math/bcl_math_histogram_2d.h"
#include "math/bcl_math_running_average_sd.h"
#include "math/bcl_math_statistics.h"
#include "math/bcl_math_template_instantiations.h"
#include "pdb/bcl_pdb_factory.h"
#include "pdb/bcl_pdb_handler.h"
#include "restraint/bcl_restraint_atom_distance.h"
#include "restraint/bcl_restraint_handler_atom_distance_assigned.h"
#include "score/bcl_score_restraint_distance_spin_label.h"

namespace bcl
{
  namespace app
  {
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class StatisticSpinLabel
    //! @brief app is for doing statistics involving spin label approximations.
    //!
    //! @author alexanns
    //! @date 09/08/08
    //!
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API StatisticSpinLabel :
      public Interface
    {

    private:

    //////////
    // data //
    //////////

      //! flag for limiting the residues that are considered so they aren't too buried
      util::ShPtr< command::FlagStatic> m_LimitNeighbors;
      util::ShPtr< command::ParameterInterface> m_NeighborLimitMeasure;         //!<param for specifying the neighbor measure to be used
      util::ShPtr< command::ParameterInterface> m_NeighborLimit;                //!< threshold value for determining if a resi will be used
      util::ShPtr< command::ParameterInterface> m_NeighborLimitSeqExcl;         //!< amount of seq. excl. that should be used
      util::ShPtr< command::ParameterInterface> m_CalculateOnlyAccessibilities; //!< only calculate accessibility values

      //! flag for specifying the output file path and name
      util::ShPtr< command::FlagInterface> m_HistogramOutputFile;

      //! flag for specifying the histogram which contains the statistics
      util::ShPtr< command::FlagStatic> m_HistogramSpecifications;
      util::ShPtr< command::ParameterInterface> m_LowerLimit;   //!< lowerlimit of the histogram
      util::ShPtr< command::ParameterInterface> m_BinSize;      //!< size of the bins of the histogram
      util::ShPtr< command::ParameterInterface> m_NumberOfBins; //!< the number of bins the histogram should contain

      //! flag for specifying spin label length
      util::ShPtr< command::FlagStatic> m_SpinLabelLength;

      //! flag for specifying the minimum CA->CB->SL angle
      util::ShPtr< command::FlagStatic> m_MinCACBSLAngle;

      //! flag for specifying the maximum CA->CB->SL angle
      util::ShPtr< command::FlagStatic> m_MaxCACBSLAngle;

      //! flag for specifying the minimum dihedral angle of the spin label
      util::ShPtr< command::FlagStatic> m_MinSLDihedralAngle;

      //! flag for specifying the maximum dihedral angle of the spin label
      util::ShPtr< command::FlagStatic> m_MaxSLDihedralAngle;

      //! flag for just calculating the radius of gyration of a protein model and nothing else
      util::ShPtr< command::FlagDynamic> m_CalculateRadiusOfGyration;

      //! flag for specifying that (Dsl-Dcb) - (Dsl-Dcb) statistics should be calculated
      util::ShPtr< command::FlagStatic> m_CalculateSL_CB_SL_CB;
      util::ShPtr< command::ParameterInterface> m_StartOffset;          //!< how far away from initial position should calcs start
      util::ShPtr< command::ParameterInterface> m_NumberOFSteps;        //!< how many steps should be made, starting at "m_Offset"
      util::ShPtr< command::ParameterInterface> m_CaCbSlAngleDrift;     //!<how much should walking SL be able to differ from initial
      util::ShPtr< command::ParameterInterface> m_SlDihedralAngleDrift; //!< amount walking SL be able to differ from initial
      util::ShPtr< command::ParameterInterface> m_RecordDslVsDCB;       //!<store Dslinitial-Dslwalking & Dcbinitial-Dcbwalking values?

      //! flag for specifying that (Dsl-Dcb) statistics should be calculated
      util::ShPtr< command::FlagStatic> m_CalculateSL_CB;
      util::ShPtr< command::ParameterInterface> m_LowerLimitX;   //!< lowerlimit of the histogram2D in the X direction
      util::ShPtr< command::ParameterInterface> m_LowerLimitY;   //!< lowerlimit of the histogram2D in the Y direction
      util::ShPtr< command::ParameterInterface> m_BinSizeX;      //!< size of the bins of the histogram2D in the X direction
      util::ShPtr< command::ParameterInterface> m_BinSizeY;      //!< size of the bins of the histogram2D in the Y direction
      util::ShPtr< command::ParameterInterface> m_NumberOfBinsX; //!< the number of bins in the histogram2D in the X direction
      util::ShPtr< command::ParameterInterface> m_NumberOfBinsY; //!< the number of bins in the histogram2D in the Y direction

      util::ShPtr< command::FlagStatic> m_CalculateExperimentalSL_CB;
      util::ShPtr< command::ParameterInterface> m_RestraintFile;

      //! flag for specifying the types of sses for Dsl-Dcb-Dsl-Dcb calculations for the spin label that is not stepped
      util::ShPtr< command::FlagInterface> m_SSETypes_Constant;

      //! flag for specifying the types of sses for Dsl-Dcb-Dsl-Dcb calculations for the spin label that is stepped
      util::ShPtr< command::FlagInterface> m_SSETypes_Step;

      //! flag for specifying that SL neighbor count minus CB neighbor count statistics should be calculated
      util::ShPtr< command::FlagStatic> m_CalculateAccessibilityStatistics;
      util::ShPtr< command::ParameterInterface> m_AccessibilityMeasure;        //!< specifies the neighbor measure that should be used
      util::ShPtr< command::ParameterInterface> m_AccessibilityMeasureSeqExcl; //!< specifies the seq. excl. that should be used

      //! flag for specifying a dihedral angle should be calculated A->B-x->C->D
      util::ShPtr< command::FlagStatic> m_CalculateDihedralAngle;
      //! flag for specifying a projection angle should be calculated B->AC->D
      util::ShPtr< command::FlagStatic> m_CalculateProjectionAngle;
      util::ShPtr< command::ParameterInterface> m_XCoordinateA;
      util::ShPtr< command::ParameterInterface> m_YCoordinateA;
      util::ShPtr< command::ParameterInterface> m_ZCoordinateA;
      util::ShPtr< command::ParameterInterface> m_XCoordinateB;
      util::ShPtr< command::ParameterInterface> m_YCoordinateB;
      util::ShPtr< command::ParameterInterface> m_ZCoordinateB;
      util::ShPtr< command::ParameterInterface> m_XCoordinateC;
      util::ShPtr< command::ParameterInterface> m_YCoordinateC;
      util::ShPtr< command::ParameterInterface> m_ZCoordinateC;
      util::ShPtr< command::ParameterInterface> m_XCoordinateD;
      util::ShPtr< command::ParameterInterface> m_YCoordinateD;
      util::ShPtr< command::ParameterInterface> m_ZCoordinateD;

      //! flag for specifying that a histogram 2d should be created from a file of value pairs
      util::ShPtr< command::FlagStatic> m_CreateHistogram2D;
      util::ShPtr< command::ParameterInterface> m_Histogram2DValuePairFilename;     //!< file containing histogram 2d values
      util::ShPtr< command::ParameterInterface> m_Histogram2DOutputFilename;        //!< file histogram 2d should be written to
      util::ShPtr< command::ParameterInterface> m_Histogram2DGnuplotOutputFilename; //!< gnuplot file of histogram 2d to write to
      util::ShPtr< command::ParameterInterface> m_Histogram2DMinX;                  //!< minimum starting value of x bins
      util::ShPtr< command::ParameterInterface> m_Histogram2DMinY;                  //!< minimum starting value of y bins
      util::ShPtr< command::ParameterInterface> m_Histogram2DBinSizeX;              //!< size of x bins
      util::ShPtr< command::ParameterInterface> m_Histogram2DBinSizeY;              //!< size of y bins
      util::ShPtr< command::ParameterInterface> m_Histogram2DNumBinsX;              //!< number of x bins
      util::ShPtr< command::ParameterInterface> m_Histogram2DNumBinsY;              //!< number of y bins
      util::ShPtr< command::ParameterInterface> m_Histogram2DGnuplotTitle;          //!< the title the gnuplot will have
      //! the column in "m_Histogram2DValuePairFilename" holding the x values
      util::ShPtr< command::ParameterInterface> m_Histogram2DValueColumnX;
      //! the column in "m_Histogram2DValuePairFilename" holding the y values
      util::ShPtr< command::ParameterInterface> m_Histogram2DValueColumnY;
      //! flag determining if histogram is normalized or not
      util::ShPtr< command::ParameterInterface> m_Histogram2DNormalize;

      //! flag for specifying that a histogram should be created from a file of value pairs
      util::ShPtr< command::FlagStatic> m_CreateHistogram;
      util::ShPtr< command::ParameterInterface> m_HistogramValueFilename;         //!< file containing histogram 2d values
      util::ShPtr< command::ParameterInterface> m_HistogramOutputFilename;        //!< file histogram 2d should be written to
      util::ShPtr< command::ParameterInterface> m_HistogramGnuplotOutputFilename; //!< gnuplot file of histogram 2d to write to
      util::ShPtr< command::ParameterInterface> m_HistogramMin;                   //!< minimum starting value of x bins
      util::ShPtr< command::ParameterInterface> m_HistogramBinSize;               //!< size of x bins
      util::ShPtr< command::ParameterInterface> m_HistogramNumBins;               //!< number of x bins
      util::ShPtr< command::ParameterInterface> m_HistogramGnuplotTitle;          //!< the title the gnuplot will have
      //! the column in "m_HistogramValueFilename" holding the x values
      util::ShPtr< command::ParameterInterface> m_HistogramValueColumn;
      //! the column in "m_HistogramValueFilename" holding the weights of x values
      util::ShPtr< command::ParameterInterface> m_HistogramWeightColumn;
      //! flag determining if histogram is normalized or not
      util::ShPtr< command::ParameterInterface> m_HistogramNormalize;

      //! flag for indicating that spin labels should be placed so that they are not clashing with other atoms
      util::ShPtr< command::FlagStatic> m_NoSLClashing;
      util::ShPtr< command::ParameterInterface> m_ClashThreshold;
      util::ShPtr< command::ParameterInterface> m_ClashMaxTrials;

      //! flag for calculating agreement with KB SpinLabel distance potential
      util::ShPtr< command::FlagStatic> m_CalculateEPRDistanceAgreement;
      util::ShPtr< command::ParameterInterface> m_DistanceRestraintFile;          //!< file containing the epr distances
      util::ShPtr< command::ParameterInterface> m_DistanceRestraintScoreOutputFile; //!< file to output scores to
      util::ShPtr< command::ParameterInterface> m_EPRDistanceAgreementPDBList;    //!< pdb file to score with restraints
      util::ShPtr< command::ParameterInterface> m_CalculateDistanceDifferences;   //!< give score or distance difference

      //! flag to calculate statistics for SL->H - CB->H
      util::ShPtr< command::FlagStatic> m_CalcSL_HStatistics;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      StatisticSpinLabel();

    public:

      //! @brief Clone function
      //! @return pointer to new Quality
      StatisticSpinLabel *Clone() const
      {
        return new StatisticSpinLabel( *this);
      }

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const
      {
        return GetStaticClassName( *this);
      }

    ////////////////
    // operations //
    ////////////////

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

      //! @brief GetAAsWithAcceptableNeighborCount gives a list of amino acids which match the neighbor count criteria
      //! @param PROTEIN_MODEL the list of all possible amino acids whose neighbor count will be checked
      //! @return returns a SiPtrVector of AABase which have acceptable neighbor countss
      const util::SiPtrVector< const biol::AABase> GetAAsWithAcceptableNeighborCount
      (
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

      //! @brief NeighborMeasureValueMeetsCriteria is for determining if a calculated neighbor measure meets a threshold
      //!        It indicates of a residue is buried enough
      //!        Important because often it is desired to only consider residues that are exposed during statistics
      //!        The neighbor measure value could be neighbor count or neighbor vector and depending on which it is
      //!        the check to see if a residue is exposed enough is inverted
      //!        i.e. (large neighbor count = buried) while (large neighbor vector = exposed)
      //! @param NEIGHBOR_MEASURE_VALUE value which will be checked for meeting the command line threshold for exposure
      //! @return bool true if "NEIGHBOR_MEASURE_VALUE" indicates the residue is exposed enough - false otherwise
      bool NeighborMeasureValueMeetsCriteria( const double &NEIGHBOR_MEASURE_VALUE) const;

      //! @brief GetNeighborMeasureValue calculates a neighbor measure value for a residue based on command line measure
      //!        The neighbor measure could be neighbor count or neighbor vector; this is determined by command line
      //! @param RESIDUE_AND_NEIGHBOR_LIST the residue of interest and its neighbor list
      //! @param PRINT_ACCESSIBILITY_STATISTICS boolean indicating if each value should be printed
      //! @return double which is the neighbor measure value
      double GetNeighborMeasureValue
      (
        const assemble::AANeighborList &RESIDUE_AND_NEIGHBOR_LIST,
        const bool PRINT_ACCESSIBILITY_STATISTICS = false
      ) const;

      //! AddSpinLabelsAndGiveDistanceBetween calculates the coordinates of two spin labels based on the CA and
      //! CB atoms of the amino acids provided and the center of the protein and a random CA->CB->SL angle and a random
      //! CENTER->CA->CB->SL dihedral angle
      //! @param AMINO_ACID_A biol::AABackBone for which a spin label is desired to be added
      //! @param AMINO_ACID_B biol::AABackBone for which a spin label is desired to be added
      //! @param PROTEIN the protein model for which the spin labels are being added
      //! @return returns the distance between the spin labels added to "AMINO_ACID_A" and "AMINO_ACID_B"
      double AddSpinLabelsAndGiveDistanceBetween
      (
        const biol::AABase &AMINO_ACID_A,
        const biol::AABase &AMINO_ACID_B,
        const assemble::ProteinModel &PROTEIN
      ) const;

      //! @brief GetSpinLabelCoordinates calculates the coordinates of s spin label based on the CA and
      //! CB atoms of the amino acid provided and the center of the protein and a random CA->CB->SL angle and a random
      //! CENTER->CA->CB->SL dihedral angle
      //! @param AMINO_ACID biol::AABackBone for which a spin label is desired to be added
      //! @param PROTEIN the protein model for which the spin labels are being added
      //! @param SL_LENGTH the length to be used for the spin label
      //! @param MIN_CA_CB_SL_ANGLE the minimum possible value for the CA->CB->SL angle
      //! @param MAX_CA_CB_SL_ANGLE the maximum possible value for the CA->CB->SL angle
      //! @param MIN_SL_DIHEDRAL_ANGLE the minimum possible value for the spin label dihedral angle
      //! @param MAX_SL_DIHEDRAL_ANGLE the minimum possible value for the spin label dihedral angle
      //! @return returns the coordinates of the spin label added to "AMINO_ACID" and the angles used to make it
      //!         The first double in the VectorND is the CA->CB->SL angle.
      //!         The second double in the VectorND is the ProtCent->CA->CB->SL dihedral angle.
      const storage::Pair< storage::VectorND< 2, double>, linal::Vector3D> GetSpinLabelCoordinates
      (
        const biol::AABase &AMINO_ACID,
        const assemble::ProteinModel &PROTEIN,
        const double &SL_LENGTH,
        const double &MIN_CA_CB_SL_ANGLE,
        const double &MAX_CA_CB_SL_ANGLE,
        const double &MIN_SL_DIHEDRAL_ANGLE,
        const double &MAX_SL_DIHEDRAL_ANGLE
      ) const;

      //! @brief CalculateRadiusOfGyration calculate the radius of gyration of a protein model
      //! @param PROTEIN_MODEL is the protein model whose radius of gyration is desired to be calculated
      //! @return return double which is the radius of gyration
      double CalculateRadiusOfGyration( const assemble::ProteinModel &PROTEIN_MODEL) const;

      //! @brief CalculateSL_CB_SL_CBStatistics does statistics for difference between two SL-CB measurements
      //! (Dsl-Dcb)-(Dsl-Dcb)
      //! @param AA_LIST_A TODO: document
      //! @param AA_LIST_B TODO: document
      //! @param NEIGHBOR_LIST TODO: document
      //! @param OFFSET_HISTOGRAM_VECTOR the Vector which contains the histograms for each step
      //! @param TOTAL_HISTOGRAM the histogram which all the statistics for AA_LIST will be inserted
      //! @param PROTEIN_MODEL is the protein model which "AA_LIST" corresponds to
      //! @param DSL_VS_DCB_POINTS TODO: document
      template< typename t_DataType_A, typename t_DataType_B>
      void CalculateSL_CB_SL_CBStatistics
      (
        const storage::Vector< t_DataType_A> &AA_LIST_A,
        const storage::Vector< t_DataType_B> &AA_LIST_B,
        const assemble::AANeighborListContainer &NEIGHBOR_LIST,
        storage::Vector< math::Histogram> &OFFSET_HISTOGRAM_VECTOR,
        math::Histogram &TOTAL_HISTOGRAM,
        const assemble::ProteinModel &PROTEIN_MODEL,
        storage::Vector< storage::VectorND< 2, storage::Vector< double> > > &DSL_VS_DCB_POINTS
      ) const;

      //! @brief CalculateSL_CBDistance takes two AABases and returns the Dsl-Dcb value for them
      //! @param AA_BASE_A is the first AABase involved in the calculation
      //! @param AA_BASE_B is the second AABase involved in the calculation
      //! @param PROTEIN_MODEL is the protein model in which "AA_BASE_A" and "AA_BASE_B" are
      //! @return double which is the Dsl-Dcb value for "AA_BASE_A" and "AA_BASE_B"
      double CalculateSL_CBDistance
      (
        const biol::AABase &AA_BASE_A, const biol::AABase &AA_BASE_B, const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

      //! @brief WriteHistogram outputs a math::Histogram to a specified file
      //! @param FILENAME the file which the math::Histogram will be output to
      //! @param HISTOGRAM is the histogram which will be output to "FILENAME"
      //! @return void
      template< typename t_HistogramType> void WriteHistogram( const std::string &FILENAME, const t_HistogramType &HISTOGRAM) const;

      //! @brief CalculateSL_CBStatistics does statistics for difference between SL and CB (Dsl-Dcb)
      //! @param AA_LIST SiPtrVector of amino acids for the current protein model that should be used
      //! @param TOTAL_HISTOGRAM the histogram which all the statistics for AA_LIST will be inserted
      //! @param HISTOGRAM_2D TODO: document
      void CalculateSL_CBStatistics
      (
        const util::SiPtrVector< const biol::AABase> &AA_LIST,
        math::Histogram &TOTAL_HISTOGRAM,
        const assemble::ProteinModel &PROTEIN_MODEL,
        math::Histogram2D &HISTOGRAM_2D
      ) const;

      //! @brief CalculateExperimentalSL_CBStatistics does statistics for difference between SL and CB (Dsl-Dcb)
      //! @param RESTRAINT_FILENAME TODO: document
      //! @param TOTAL_HISTOGRAM the histogram which all the statistics for AA_LIST will be inserted
      //! @param PROTEIN_MODEL TODO: document
      //! @param HISTOGRAM_2D TODO: document
      void CalculateExperimentalSL_CBStatistics
      (
        const std::string &RESTRAINT_FILENAME,
        math::Histogram &TOTAL_HISTOGRAM,
        const assemble::ProteinModel &PROTEIN_MODEL,
        math::Histogram2D &HISTOGRAM_2D
      ) const;

      //! @brief GetSSETypesFromCommandLine converts the SSE strings given over the command line into a Set of SSTypes
      //! @param FLAG is the command::FlagInterface which has the list of SSE strings
      //! @return Set of biol::SSTypes which were in the "FLAG"
      const storage::Set< biol::SSType> GetSSETypesFromCommandLine( const command::FlagInterface &FLAG) const;

      //! @brief SL_CB_SL_CBStatisticsSSESpecific does (SL-CB)-(SL-CB) statistics using specific SSE types
      //! @param PROTEIN_MODEL the protein model for which the (SL-CB)-(SL-CB) statistics will be done
      //! @param OFFSET_HISTOGRAM_VECTOR the Vector which contains the histograms for each step
      //! @param TOTAL_HISTOGRAM the histogram which all the statistics for AA_LIST will be inserted
      //! @param DSL_VS_DCB_POINTS TODO: document
      void SL_CB_SL_CBStatisticsSSESpecific
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        storage::Vector< math::Histogram> &OFFSET_HISTOGRAM_VECTOR,
        math::Histogram &TOTAL_HISTOGRAM,
        storage::Vector< storage::VectorND< 2, storage::Vector< double> > > &DSL_VS_DCB_POINTS
      ) const;

      //! @brief CalculateSL_CBStatistics does statistics for difference between SL and CB (Dsl-Dcb)
      //! @param AA_LIST SiPtrVector of amino acids for the current protein model that should be used
      //! @param TOTAL_HISTOGRAM the histogram which all the statistics for AA_LIST will be inserted
      //! @param PROTEIN_MODEL TODO: document
      void CalculateSLNC_CBNCStatistics
      (
        const util::SiPtrVector< const biol::AABase> &AA_LIST,
        math::Histogram &TOTAL_HISTOGRAM,
        const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

      //! @brief CalculateCBCACAAngleSummation is for calculating the sum of the CB->CA->CA angle between two resis
      //!        the angle is calculated going from CBa->CAa->CAb and from CBb->CAb->CAa and these two values summed
      //!        this can be used to see if there is a correlation between the Dsl-Dcb value and the relative
      //!        orientations of the two spin labels
      //! @param AA_BASE_A the first amino acid involved in the calculation
      //! @param AA_BASE_B the first amino acid involved in the calculation
      //! @return double indicating the radian sum of the two angles CBa->CAa->CAb and CBb->CAb->CAa
      double CalculateCBCACAAngleSummation( const biol::AABase &AA_BASE_A, const biol::AABase &AA_BASE_B) const;

      //! @brief CreateHistogram2DFromValuePairFile is for creating a histogram 2d and gnuplot script from a file
      void CreateHistogram2DFromValuePairFile() const;

      //! @brief CreateHistogramFromValueFile is for creating a histogram and gnuplot script from a file
      void CreateHistogramFromValueFile() const;

      //! @brief SLClash is for determining if a psuedo-spin label clashes with other atoms
      //! @param SL_COORDINATES the coordinates of the psuedo-spin label
      //! @param PROTEIN the protein in which the psuedo-spin label was placed
      //! @param CLASH_THRESHOLD the distance between the psuedo-spin label and other atoms which qualifies as a clash
      //! @return bool true if there is a clash - false if the psuedo spin label does not clash with other atoms
      bool SLClash
      (
        const linal::Vector3D &SL_COORDINATES, const assemble::ProteinModel &PROTEIN, const double CLASH_THRESHOLD
      ) const;

      //! @brief CalculateEPRDistanceAgreement for calculating the agreement of a protein model with distance restraints
      void CalculateEPRDistanceAgreement() const;

      //! @brief GetProteinModel is for creating a protein model given a pdb filename
      //! @param PDB_FILENAME is the name of the pdb filename from which a protein model will be created
      //! @return assemble::ProteinModel created from "PDB_FILENAME"
      assemble::ProteinModel GetProteinModel( const std::string &PDB_FILENAME) const;

      //! @brief function to calculate a value for SL->H - CB->H
      //! @param AA_BASE_A the residue whose amide hydrogen will be used in the current statistic
      //! @param AA_BASE_B the residue whose cb will be used in the current statistic and will have pseudo SL added
      //! @param PROTEIN_MODEL mode where the residues are located
      //! @return double which gives the value for SL->H - CB->H
      double CalculateSL_NHDistance
      (
        const biol::AABase &AA_BASE_A, const biol::AABase &AA_BASE_B, const assemble::ProteinModel &PROTEIN_MODEL
      ) const;

    ///////////////
    // operators //
    ///////////////

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief read from std::istream
      //! @param ISTREAM input stream
      //! @return istream which was read from
      std::istream &Read( std::istream &ISTREAM)
      {
        return ISTREAM;
      }

      //! @brief write to std::ostream
      //! @param OSTREAM output stream to write to
      //! @param INDENT number of indentations
      //! @return output stream which was written to
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const
      {
        return OSTREAM;
      }

    //////////////////////
    // helper functions //
    //////////////////////

    private:

      static const ApplicationType StatisticSpinLabel_Instance;

    }; // StatisticSpinLabel

    //! @brief initializes the command object for that executable
    util::ShPtr< command::Command> StatisticSpinLabel::InitializeCommand() const
    {
      util::ShPtr< command::Command> sp_cmd( new command::Command());

      sp_cmd->AddFlag( assemble::ProteinEnsemble::GetFlagPDBList());

      // flag for specifying the neighbor limit for residues
      sp_cmd->AddFlag( m_LimitNeighbors);

      // input m_HistogramOutputFile
      sp_cmd->AddFlag( m_HistogramOutputFile);

      // input m_HistogramSpecifications
      sp_cmd->AddFlag( m_HistogramSpecifications);

      // flag for specifying spin label length
      sp_cmd->AddFlag( m_SpinLabelLength);

      // flag for specifying the minimum CA->CB->SL angle
      sp_cmd->AddFlag( m_MinCACBSLAngle);

      // flag for specifying the maximum CA->CB->SL angle
      sp_cmd->AddFlag( m_MaxCACBSLAngle);

      // flag for specifying the minimum dihedral angle of the spin label
      sp_cmd->AddFlag( m_MinSLDihedralAngle);

      // flag for specifying the maximum dihedral angle of the spin label
      sp_cmd->AddFlag( m_MaxSLDihedralAngle);

      // add flag for specifying the minimum lengths of SSEs to be included in the protein models
      sp_cmd->AddFlag( pdb::Factory::GetFlagMinSSESize());

      // add flag for calculating the radius of gyration of a protein model
      sp_cmd->AddFlag( m_CalculateRadiusOfGyration);

      // add flag for specifying that (Dsl-Dcb) - (Dsl-Dcb) statistics should be calculated
      sp_cmd->AddFlag( m_CalculateSL_CB_SL_CB);

      // add flag for specifying that (Dsl-Dcb) statistics should be calculated
      sp_cmd->AddFlag( m_CalculateSL_CB);

      sp_cmd->AddFlag( m_CalculateExperimentalSL_CB);

      // flag for specifying the types of sses for Dsl-Dcb-Dsl-Dcb calculations for the spin label that is not stepped
      sp_cmd->AddFlag( m_SSETypes_Constant);

      // flag for specifying the types of sses for Dsl-Dcb-Dsl-Dcb calculations for the spin label that is not stepped
      sp_cmd->AddFlag( m_SSETypes_Step);

      // flag for doing accessibility statistics
      sp_cmd->AddFlag( m_CalculateAccessibilityStatistics);

      // flag for calculating a dihedral angle
      sp_cmd->AddFlag( m_CalculateDihedralAngle);

      // flag for calculating a projection angle
      sp_cmd->AddFlag( m_CalculateProjectionAngle);

      // flag for specifying that a histogram 2d should be created from a file of value pairs
      sp_cmd->AddFlag( m_CreateHistogram2D);

      // flag for specifying that a histogram should be created from a file of values
      sp_cmd->AddFlag( m_CreateHistogram);

      sp_cmd->AddFlag( m_NoSLClashing);

      // flag for specifying that agreement with restraints should be calculated
      sp_cmd->AddFlag( m_CalculateEPRDistanceAgreement);

      // flag to calculate statistics for SL->H - CB->H
      sp_cmd->AddFlag( m_CalcSL_HStatistics);

      // add flag to specify the amino acid class that should be used when creating the protein model
      sp_cmd->AddFlag( pdb::Factory::GetFlagAAClass());

      // add default bcl parameters
      command::GetAppDefaultFlags().AddDefaultCommandlineFlags( *sp_cmd);

      // return assembled Command object
      return sp_cmd;
    }

    //! @brief the Main function
    //! @return error code - 0 for success
    int StatisticSpinLabel::Main() const
    {
      // true if the "m_CalculateDihedralAngle" flag was set
      if( m_CalculateDihedralAngle->GetFlag() || m_CalculateProjectionAngle->GetFlag())
      {
        // create Vector3d "point_a" and initialize with the first set of coordinates
        const linal::Vector3D point_a
        (
          m_XCoordinateA->GetNumericalValue< double>(),
          m_YCoordinateA->GetNumericalValue< double>(),
          m_ZCoordinateA->GetNumericalValue< double>()
        );

        // create Vector3d "point_b" and initialize with the second set of coordinates
        const linal::Vector3D point_b
        (
          m_XCoordinateB->GetNumericalValue< double>(),
          m_YCoordinateB->GetNumericalValue< double>(),
          m_ZCoordinateB->GetNumericalValue< double>()
        );

        // create Vector3d "point_c" and initialize with the third set of coordinates
        const linal::Vector3D point_c
        (
          m_XCoordinateC->GetNumericalValue< double>(),
          m_YCoordinateC->GetNumericalValue< double>(),
          m_ZCoordinateC->GetNumericalValue< double>()
        );

        // create Vector3d "point_d" and initialize with the fourth set of coordinates
        const linal::Vector3D point_d
        (
          m_XCoordinateD->GetNumericalValue< double>(),
          m_YCoordinateD->GetNumericalValue< double>(),
          m_ZCoordinateD->GetNumericalValue< double>()
        );

        // create ofstream "write" for writing radius of gyration to a file
        io::OFStream write;

        if( m_CalculateDihedralAngle->GetFlag())
        {
          // calculate the dihedral angle
          const double dihedral_angle( linal::Dihedral( point_a, point_b, point_c, point_d));

          // open "write" and bind it to "m_HistogramOutputFile"
          io::File::MustOpenOFStream
          (
            write, m_HistogramOutputFile->GetFirstParameter()->GetValue() + "dihedral"
          );

          // output "dihedral_angle" to "write"
          write << dihedral_angle << '\n';
          io::File::CloseClearFStream( write);
        }

        if( m_CalculateProjectionAngle->GetFlag())
        {
          // create double "projection_angle" and initialize with the projection angle between the four points
          // angle is between b->ac->d
          const double projection_angle( linal::ProjAngle( point_a, point_b, point_c, point_d));
          // open "write" and bind it to "m_HistogramOutputFile"
          io::File::MustOpenOFStream
          (
            write, m_HistogramOutputFile->GetFirstParameter()->GetValue() + "projection_angle"
          );

          write << projection_angle;
          io::File::CloseClearFStream( write);
        }

        return 0;
      }

      // true if "m_CreateHistogram2D" flag was set in command line - just create histogram 2d from data then exit
      if( m_CreateHistogram2D->GetFlag())
      {
        CreateHistogram2DFromValuePairFile();
        return 0;
      }

      // true if "m_CreateHistogram" flag was set in command line - just create histogram from data then exit
      if( m_CreateHistogram->GetFlag())
      {
        CreateHistogramFromValueFile();
        return 0;
      }

      // make sure that if the "m_CalculateSL_CB_SL_CB" flag was set, then also the desired SSEs were provided also
      if( m_CalculateSL_CB_SL_CB->GetFlag())
      {
        BCL_Assert
        (
          m_SSETypes_Constant->GetFlag() && m_SSETypes_Step->GetFlag(),
          "When using the flag calculateSL_CB_SL_CB, you must explicitly specify they types of SSEs to be used for constant and walking positions by using the \"sse_types_step\" and \"sse_types_constant flags\"."
        );
      }
      // make sure that if the SSETypes flags were set, then also the "m_CalculateSL_CB_SL_CB" flag was set
      if( m_SSETypes_Constant->GetFlag() || m_SSETypes_Step->GetFlag())
      {
        BCL_Assert
        (
          m_CalculateSL_CB_SL_CB->GetFlag(),
          "Using the flags \"sse_types_step\" and/or \"sse_types_constant\" only work with the \"calculateSL_CB_SL_CB\" flag."
        );
      }

      //create Histogram "histogram" for holding statistics and initialize with command line parameters
      math::Histogram histogram
      (
        m_LowerLimit->GetNumericalValue< double>(),
        m_BinSize->GetNumericalValue< double>(),
        m_NumberOfBins->GetNumericalValue< size_t>()
      );

      math::Histogram2D histogram_2d
      (
        storage::VectorND< 2, double>
        (
          m_LowerLimitX->GetNumericalValue< double>(),
          m_LowerLimitY->GetNumericalValue< double>()
        ),
        storage::VectorND< 2, double>
        (
          m_BinSizeX->GetNumericalValue< double>(),
          m_BinSizeY->GetNumericalValue< double>()
        ),
        storage::VectorND< 2, size_t>
        (
          m_NumberOfBinsX->GetNumericalValue< size_t>(),
          m_NumberOfBinsY->GetNumericalValue< size_t>()
        )
      );

      // create storage::Vector "histogram_vector" and initialize with a "histogram" for each step
      storage::Vector< math::Histogram> histogram_vector
      (
        m_NumberOFSteps->GetNumericalValue< int>(), histogram
      );

      // create math::Histogram "sl_cb_sl_cb_histogram"for holding all the statistics related to Dsl-Dcb - Dsl-Dcb
      // calculations and initialize with "histogram"
      math::Histogram sl_cb_sl_cb_histogram( histogram);

      // create Vector "dsl_vs_dcb_points" for holding data if the "m_RecordDslVsDCB" parameter is true
      storage::Vector< storage::VectorND< 2, storage::Vector< double> > > dsl_vs_dcb_points
      (
        m_NumberOFSteps->GetNumericalValue< int>(),
        storage::VectorND< 2, storage::Vector< double> >( storage::Vector< double>(), storage::Vector< double>())
      );

      // true if only agreement with experimental epr distances is desired to be calculated
      if( m_CalculateEPRDistanceAgreement->GetFlag())
      {
        CalculateEPRDistanceAgreement();
        return 0;
      }

      // create Vector "pdb_list" initialize with the list of strings obtained from contained in the pdb list file
      const storage::Vector< std::string> pdb_list
      (
        command::StringVectorFromFilenameParameter( *assemble::ProteinEnsemble::GetFlagPDBList()->GetFirstParameter())
      );

      // loop through "file_path_and_name" to do statistics over all proteins
      for
      (
        storage::Vector< std::string>::const_iterator file_itr( pdb_list.Begin()), file_itr_end( pdb_list.End());
        file_itr != file_itr_end; ++file_itr
      )
      {
        // create string "pdb_filename" and initialize with filename currently denoted by "itr"
        const std::string pdb_filename( *file_itr);

        // create ProteinModel "protein_model" and initialize with protein model created from "pdb_filename"
        const assemble::ProteinModel protein_model( GetProteinModel( pdb_filename));

        // message that the protein model "pdb_filename" was created
        BCL_MessageCrt( "protein model " + pdb_filename + " created");

        // true if it is desired to only calculate the radius of gyration
        if( m_CalculateRadiusOfGyration->GetFlag())
        {
          // create ofstream "write" for writing radius of gyration to a file
          io::OFStream write;

          // open "write" and bind it to "m_HistogramOutputFile"
          io::File::MustOpenOFStream
          (
            write, m_HistogramOutputFile->GetFirstParameter()->GetValue() + pdb_filename + "rgyr"
          );

          write << pdb_filename << " radius_gyration " << CalculateRadiusOfGyration( protein_model) << '\n';
          io::File::CloseClearFStream( write);

          return 0;
        }

        // create ShPtrVector "all_chains" and initialize with the chains in "protein_model"
        const util::ShPtrVector< assemble::Chain> all_chains( protein_model.GetChains());

        // create SiPtrVector "all_amino_acids" and initialize with the amino acids in the SSEs of "protein_model"
        const util::SiPtrVector< const biol::AABase> all_amino_acids( protein_model.GetAminoAcids());

        // print the number of amino acids in "all_amino_acids"
        BCL_MessageStd
        (
          "The number of amino acids in " + pdb_filename + " is " + util::Format()( all_amino_acids.GetSize())
        );

        // create SiPtrVector "exposed_amino_acids" and initialize with the amino acids which have an acceptable
         // neighbor count
         util::SiPtrVector< const biol::AABase> exposed_amino_acids;

        if( m_LimitNeighbors->GetFlag())
        {
          // create SiPtrVector "exposed_amino_acids" and initialize with the amino acids which have an acceptable
          // neighbor count
          exposed_amino_acids = GetAAsWithAcceptableNeighborCount( protein_model);

          if( m_CalculateOnlyAccessibilities->GetValue() == "true")
          {
            continue;
          }
        }
        else
        {
          exposed_amino_acids = all_amino_acids;
        }

        // true if "m_CalculateSL_CB" was set in the command line
        if( m_CalculateSL_CB->GetFlag() || m_CalcSL_HStatistics->GetFlag())
        {
          BCL_MessageCrt( "Calculating only CalculateSL_CBStatistics");
          if( m_CalculateExperimentalSL_CB->GetFlag())
          {
            BCL_MessageCrt( "Calculating only ExperimentalCalculateSL_CBStatistics");
            // calculate ExperimentalDsl-Dcb statistics
            CalculateExperimentalSL_CBStatistics
            (
              m_RestraintFile->GetValue(),
              histogram,
              protein_model,
              histogram_2d
            );
          }
          else
          {
            // calculate Dsl-Dcb statistics
            CalculateSL_CBStatistics( exposed_amino_acids, histogram, protein_model, histogram_2d);
          }
        }
        else if( m_CalculateAccessibilityStatistics->GetFlag())
        {
          BCL_MessageCrt( "Calculating only AccessibilityStatistics");
          CalculateSLNC_CBNCStatistics( exposed_amino_acids, histogram, protein_model);
        }
        // true if "m_CalculateSL_CB_SL_CB" was set in the command line
        else if( m_CalculateSL_CB_SL_CB->GetFlag())
        {
          BCL_MessageCrt( "Calculating only CalculateSL_CB_SL_CBStatistics");
          // calculate Dsl-Dcb - Dsl-Dcb statistics
          SL_CB_SL_CBStatisticsSSESpecific( protein_model, histogram_vector, sl_cb_sl_cb_histogram, dsl_vs_dcb_points);
        }
      }

      // true if "m_CalcSL_HStatistics" was set in the command line
      if( m_CalcSL_HStatistics->GetFlag())
      {
        // write "histogram" which contains the Dsl-Dcb statistics
        WriteHistogram( m_HistogramOutputFile->GetFirstParameter()->GetValue() + "_pre", histogram);
      }

      // true if "m_CalculateSL_CB" was set in the command line
      if( m_CalculateSL_CB->GetFlag())
      {
        // write "histogram" which contains the Dsl-Dcb statistics
        WriteHistogram( m_HistogramOutputFile->GetFirstParameter()->GetValue() + "_Dsl-Dcb", histogram);
        // write "histogram_2d" which contains the Dsl-Dcb and angle 2D statistics
        WriteHistogram( m_HistogramOutputFile->GetFirstParameter()->GetValue() + "_Dsl-DcbVSangle_2d", histogram_2d);
      }
      else if //< true if Accessibility statistics were calculated using neighbor count
      (
        m_CalculateAccessibilityStatistics->GetFlag() &&
        m_AccessibilityMeasure->GetValue() == "NeighborCount"
      )
      {
        WriteHistogram( m_HistogramOutputFile->GetFirstParameter()->GetValue() + "_Dnc-Dnc", histogram);
      }
      else if //< true if Accessibility statistics were calculated using neighbor vector
      (
        m_CalculateAccessibilityStatistics->GetFlag() &&
        m_AccessibilityMeasure->GetValue() == "NeighborVector"
      )
      {
        WriteHistogram( m_HistogramOutputFile->GetFirstParameter()->GetValue() + "_Dnv-Dnv", histogram);
      }
      // true if "m_CalculateSL_CB_SL_CB" was set in the command line
      else if( m_CalculateSL_CB_SL_CB->GetFlag())
      {
        // iterate through "histogram_vector" and write out each of the histograms
        for
        (
          storage::Vector< math::Histogram>::const_iterator
            itr( histogram_vector.Begin()), itr_end( histogram_vector.End());
          itr != itr_end;
          ++itr
        )
        {
          const size_t histogram_number( itr - histogram_vector.Begin());
          std::string hist_num_string;
          if( histogram_number < 10)
          {
            hist_num_string = "0" + util::Format()( histogram_number);
          }
          else
          {
            hist_num_string = util::Format()( histogram_number);
          }

          // create std::string "histogram_filename" and initialize with filename for current histogram
          const std::string histogram_filename
          (
            m_HistogramOutputFile->GetFirstParameter()->GetValue()
            + "_Dsl-Dcb-Dsl-Dcb" + hist_num_string
          );

          // write out the histogram denoted by "itr" to "histogram_filename"
          WriteHistogram
          (
            histogram_filename,
            *itr
          );
        }

        // write out "sl_cb_sl_cb_histogram" histogram
        WriteHistogram
        (
          m_HistogramOutputFile->GetFirstParameter()->GetValue() + "_Dsl-Dcb-Dsl-DcbTotal",
          sl_cb_sl_cb_histogram
        );

        // true if dslinitial-dslwalking and dcbinitial-dcbwalking values were stored
        if( m_RecordDslVsDCB->GetValue() == "true")
        {
          io::OFStream ofstream;
          std::string dsl_vs_dcb_filename( m_HistogramOutputFile->GetFirstParameter()->GetValue() + "_dsl_vs_dcb_correlations.txt");
          io::File::MustOpenOFStream( ofstream, dsl_vs_dcb_filename);
          // iterate through "histogram_vector" and write out each of the histograms
          for
          (
            storage::Vector< storage::VectorND< 2, storage::Vector< double> > >::const_iterator
              itr( dsl_vs_dcb_points.Begin()), itr_end( dsl_vs_dcb_points.End());
            itr != itr_end;
            ++itr
          )
          {
            const double rsquared
            (
              math::Statistics::RSquared
              (
                itr->First().Begin(),
                itr->First().End(),
                itr->Second().Begin(),
                itr->Second().End()
              )
            );
            size_t position( itr - dsl_vs_dcb_points.Begin());
            BCL_MessageDbg
            (
              "size of vector for position " + util::Format()( position) +
              " is " + util::Format()( itr->First().GetSize())
            );
            ofstream << "position : " << position << " Rsqrd : " << rsquared << '\n';
          }
          io::File::CloseClearFStream( ofstream);
        }
      }

//      GetNeighborMeasureValue
//      (
//        std::pair
//        <
//          util::SiPtr< const biol::AABase>, util::SiPtr< const assemble::AANeighborList::SingleAANeighborList>
//        >(),
//        true
//      );

      //successful end
      return 0;
    }

    //! @brief WriteHistogram outputs a math::Histogram to a specified file
    //! @param FILENAME the file which the math::Histogram will be output to
    //! @param HISTOGRAM is the histogram which will be output to "FILENAME"
    //! @return void
    template< typename t_HistogramType> void StatisticSpinLabel::WriteHistogram( const std::string &FILENAME, const t_HistogramType &HISTOGRAM) const
    {
      // print "histogram"
      BCL_MessageStd( "Histogram " + FILENAME + " is \n" + util::Format()( HISTOGRAM) + "\n");

      // create ofstream "write" for writing "histogram" to a file
      io::OFStream write;

      // open "write" and bind it to "m_HistogramOutputFile"
      io::File::MustOpenOFStream( write, FILENAME);

      // output "HISTOGRAM" to file
      write << HISTOGRAM;

      // close and clear "read"
      io::File::CloseClearFStream( write);
    }

    //! @brief CalculateSL_CBStatistics does statistics for difference between SL and CB (Dsl-Dcb)
    //! @param AA_LIST SiPtrVector of amino acids for the current protein model that should be used
    //! @param TOTAL_HISTOGRAM the histogram which all the statistics for AA_LIST will be inserted
    //! @param PROTEIN_MODEL TODO: document
    //! @param HISTOGRAM_2D TODO: document
    void StatisticSpinLabel::CalculateSL_CBStatistics
    (
      const util::SiPtrVector< const biol::AABase> &AA_LIST,
      math::Histogram &TOTAL_HISTOGRAM,
      const assemble::ProteinModel &PROTEIN_MODEL,
      math::Histogram2D &HISTOGRAM_2D
    ) const
    {
      // create dataset to calculate mean and sd
      math::RunningAverageSD< double> mean_sd;

      // loop over amino acids in "exposed_amino_acids" to get the Dsl-Dcb distances
      for
      (
        util::SiPtrVector< const biol::AABase>::const_iterator
          aa_itr_a( AA_LIST.Begin()), aa_itr_end( AA_LIST.End());
        aa_itr_a != aa_itr_end;
        ++aa_itr_a
      )
      {
        if( !( *aa_itr_a)->GetFirstSidechainAtom().AllCoordinatesDefined() || !( *aa_itr_a)->GetCA().AllCoordinatesDefined())
        {
          continue;
        }

        // loop over amino acids in "exposed_amino_acids" to get the Dsl-Dcb distances
        for
        (
          util::SiPtrVector< const biol::AABase>::const_iterator aa_itr_b( aa_itr_a + 1);
          aa_itr_b != aa_itr_end;
          ++aa_itr_b
        )
        {
          if( !( *aa_itr_b)->GetFirstSidechainAtom().AllCoordinatesDefined() || !( *aa_itr_b)->GetCA().AllCoordinatesDefined())
          {
            continue;
          }

          // true if the flag to calculate statistics between the spin label and the amide hydrogren atom was given
          if( m_CalcSL_HStatistics->GetFlag())
          {
            // calculate SL->H - H->CB statistic
            const double difference_a( CalculateSL_NHDistance( **aa_itr_a, **aa_itr_b, PROTEIN_MODEL));
            // true if difference_a is defined, need to push it back into "TOTAL_HISTOGRAM"
            if( util::IsDefined( difference_a))
            {
              TOTAL_HISTOGRAM.PushBack( difference_a);
              mean_sd += difference_a;
            }

            // calculate SL->H - H->CB statistic
            const double difference_b( CalculateSL_NHDistance( **aa_itr_b, **aa_itr_a, PROTEIN_MODEL));
            // true if difference_b is defined, need to push it back into "TOTAL_HISTOGRAM"
            if( util::IsDefined( difference_b))
            {
              TOTAL_HISTOGRAM.PushBack( difference_b);
              mean_sd += difference_b;
            }
          }
          else //< otherwise calculate the SL->SL - CB->CB statistic
          {
            const double sl_cb_difference( CalculateSL_CBDistance( **aa_itr_a, **aa_itr_b, PROTEIN_MODEL));
            const double angle_summation( CalculateCBCACAAngleSummation( **aa_itr_a, **aa_itr_b));

            // PushBack the SL-CB difference between "aa_itr_a" and "aa_itr_b" into "histogram"
            TOTAL_HISTOGRAM.PushBack( sl_cb_difference);
            HISTOGRAM_2D.PushBack( storage::VectorND< 2, double>( sl_cb_difference, angle_summation));
          }
        }
      }

      // if the flag to calculate statistics between the spin label and the amide hydrogren atom was given
      if( m_CalcSL_HStatistics->GetFlag())
      {
        BCL_MessageStd( "Counts: " + util::Format()( mean_sd.GetWeight()));
        BCL_MessageStd( "Mean: " + util::Format()( mean_sd.GetAverage()));
        BCL_MessageStd( "SD: " + util::Format()( mean_sd.GetStandardDeviation()));
      }
    }

    //! @brief CalculateExperimentalSL_CBStatistics does statistics for difference between SL and CB (Dsl-Dcb)
    //! @param TOTAL_HISTOGRAM the histogram which all the statistics for AA_LIST will be inserted
    //! @return void
    void StatisticSpinLabel::CalculateExperimentalSL_CBStatistics
    (
      const std::string &RESTRAINT_FILENAME,
      math::Histogram &TOTAL_HISTOGRAM,
      const assemble::ProteinModel &PROTEIN_MODEL,
      math::Histogram2D &HISTOGRAM_2D
    ) const
    {
      storage::Vector< storage::VectorND< 2, char> > chain_ids;
      storage::Vector< storage::VectorND< 2, size_t> > aa_ids;
      storage::Vector< storage::VectorND< 2, biol::AtomType> > atom_types;
      storage::Vector< restraint::Distance> distances;
      io::IFStream read;
      io::File::MustOpenIFStream( read, RESTRAINT_FILENAME);

      // get restraints from file
      const util::ShPtrVector< restraint::AtomDistance> restraints
      (
        restraint::HandlerAtomDistanceAssigned().ReadRestraints( read)
      );

      for
      (
        util::ShPtrVector< restraint::AtomDistance>::const_iterator
          itr( restraints.Begin()), itr_end( restraints.End());
        itr != itr_end; ++itr
      )
      {
        const util::SiPtr< const biol::AABase> aa_a( ( *itr)->GetData().First()->LocateAA( PROTEIN_MODEL));
        const util::SiPtr< const biol::AABase> aa_b( ( *itr)->GetData().Second()->LocateAA( PROTEIN_MODEL));
        const double angle_summation( CalculateCBCACAAngleSummation( *aa_a, *aa_b));
        const double distance( linal::Distance( aa_a->GetFirstSidechainAtom().GetCoordinates(), aa_b->GetFirstSidechainAtom().GetCoordinates()));

        const double dsl_dcb
        (
          ( *itr)->GetDistance()->GetDistance() - distance
        );
        TOTAL_HISTOGRAM.PushBack( dsl_dcb);
        HISTOGRAM_2D.PushBack( storage::VectorND< 2, double>( dsl_dcb, angle_summation));
        BCL_MessageDbg
        (
          "for " + ( *itr)->GetIdentification() +
          " dsl : " + util::Format()( ( *itr)->GetDistance()->GetDistance()) + " Dcb : " + util::Format()( distance)
          + " Dsl-Dcb : " + util::Format()( dsl_dcb)
        );
      }
    }

    //! @brief CalculateCBCACAAngleSummation is for calculating the sum of the CB->CA->CA angle between two resis
    //!        the angle is calculated going from CBa->CAa->CAb and from CBb->CAb->CAa and these two values summed
    //!        this can be used to see if there is a correlation between the Dsl-Dcb value and the relative
    //!        orientations of the two spin labels
    //! @param AA_BASE_A the first amino acid involved in the calculation
    //! @param AA_BASE_B the first amino acid involved in the calculation
    //! @return double indicating the radian sum of the two angles CBa->CAa->CAb and CBb->CAb->CAa
    double StatisticSpinLabel::CalculateCBCACAAngleSummation
    (
      const biol::AABase &AA_BASE_A, const biol::AABase &AA_BASE_B
    ) const
    {
      // get the angle from CBa->CAa->CAb
      const double cb_ca_ca_angle_a
      (
        linal::ProjAngle
        (
          AA_BASE_A.GetCA().GetCoordinates(),
          AA_BASE_A.GetFirstSidechainAtom().GetCoordinates(),
          AA_BASE_A.GetCA().GetCoordinates(),
          AA_BASE_B.GetFirstSidechainAtom().GetCoordinates()
        )
      );

      // get the angle from CBb->CAb->CAa
      const double cb_ca_ca_angle_b
      (
        linal::ProjAngle
        (
          AA_BASE_B.GetCA().GetCoordinates(),
          AA_BASE_B.GetFirstSidechainAtom().GetCoordinates(),
          AA_BASE_B.GetCA().GetCoordinates(),
          AA_BASE_A.GetCA().GetCoordinates()
        )
      );

      // create const double "angle_summation" and initialize with the sum of "cb_ca_ca_angle_a" and "cb_ca_ca_angle_b"
      const double angle_summation
      (
        cb_ca_ca_angle_a + cb_ca_ca_angle_b
      );

      // return "angle_summation"
      return angle_summation;
    }

    //! @brief CalculateSL_CBStatistics does statistics for difference between SL and CB (Dsl-Dcb)
    //! @param AA_LIST SiPtrVector of amino acids for the current protein model that should be used
    //! @param TOTAL_HISTOGRAM the histogram which all the statistics for AA_LIST will be inserted
    //! @return void
    void StatisticSpinLabel::CalculateSLNC_CBNCStatistics
    (
      const util::SiPtrVector< const biol::AABase> &AA_LIST,
      math::Histogram &TOTAL_HISTOGRAM,
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      static const double s_max_distance( 1000.0);
      // create AllAANeighborList "neighbor_list" initialize with the neighbors for each amino acid
      // in "AA_LIST"
      const assemble::AANeighborListContainer neighbor_list
      (
        AA_LIST, s_max_distance, m_AccessibilityMeasureSeqExcl->GetNumericalValue< int>(), true
      );

      // map to hold coords and the neighboring coords
      storage::Map< util::SiPtr< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >
        cb_coords_and_neigh_coords;

      for
      (
        assemble::AANeighborListContainer::const_iterator
          all_itr( neighbor_list.Begin()), all_itr_end( neighbor_list.End());
        all_itr != all_itr_end;
        ++all_itr
      )
      {
        const util::SiPtr< const linal::Vector3D> resi_cb_coords( all_itr->second.GetCenterAminoAcid()->GetFirstSidechainAtom().GetCoordinates());
        cb_coords_and_neigh_coords.Insert
        (
          std::pair< util::SiPtr< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D> >
          (
            resi_cb_coords, util::SiPtrVector< const linal::Vector3D>()
          )
        );
        // iterate over all of the neighbors of the AA currently denoted by "all_itr"
        for
        (
          assemble::AANeighborList::const_iterator
            neighbor_itr( all_itr->second.Begin()), neighbor_itr_end( all_itr->second.End());
          neighbor_itr != neighbor_itr_end;
          ++neighbor_itr
        )
        {
          // add the CB coords of the neighbor currently denoted by "neighbor_itr"
          cb_coords_and_neigh_coords[ resi_cb_coords].PushBack
          (
            util::ToSiPtr( neighbor_itr->First()->GetFirstSidechainAtom().GetCoordinates())
          );
        }
      }

      // create iterator "aa_itr" and iterator end "aa_itr" on "AA_LIST"
      util::SiPtrVector< const biol::AABase>::const_iterator aa_itr( AA_LIST.Begin());

      // loop over "cb_coords_and_neigh_coords" to get the SLnc-CBnc statistics
      for
      (
        storage::Map
        <
          util::SiPtr< const linal::Vector3D>, util::SiPtrVector< const linal::Vector3D>
        >::const_iterator itr( cb_coords_and_neigh_coords.Begin()), itr_end( cb_coords_and_neigh_coords.End());
        itr != itr_end;
        ++itr, ++aa_itr
      )
      {

        // true if the CA and CB coordinates are defined
        if( ( *aa_itr)->GetFirstSidechainAtom().AllCoordinatesDefined() && ( *aa_itr)->GetCA().AllCoordinatesDefined())
        {
          // create Vector3D "spin_label_coords" and initialize with the coordinates of a spin label for the residue
          // currently denoted by "aa_itr"
          const linal::Vector3D spin_label_coords
          (
            GetSpinLabelCoordinates
            (
              **aa_itr,
              PROTEIN_MODEL,
              m_SpinLabelLength->GetFirstParameter()->GetNumericalValue< double>(),
              m_MinCACBSLAngle->GetFirstParameter()->GetNumericalValue< double>(),
              m_MaxCACBSLAngle->GetFirstParameter()->GetNumericalValue< double>(),
              m_MinSLDihedralAngle->GetFirstParameter()->GetNumericalValue< double>(),
              m_MaxSLDihedralAngle->GetFirstParameter()->GetNumericalValue< double>()
            ).Second()
          );

          // true if the "m_CalculateSLNC_CBNC" flag was set in the commandline - will use neighbor count
          if( m_AccessibilityMeasure->GetValue() == "NeighborCount")
          {
            // create double "cb_nc" and initialize with calculated neighbor count for the coordinates denoted
            // by "coord_itr"
            const double cb_nc
            (
              coord::CountNeighbors
              (
                itr->second,
                *itr->first,
                assemble::AANeighborCount::GetDefaultThresholdLowHigh().GetMin(),
                assemble::AANeighborCount::GetDefaultThresholdLowHigh().GetMax()
              )
            );

            // create double "sl_nc" and initialize with calculated the neighbor count for the "spin_label_coords"
            const double sl_nc
            (
              coord::CountNeighbors
              (
                itr->second,
                spin_label_coords,
                assemble::AANeighborCount::GetDefaultThresholdLowHigh().GetMin(),
                assemble::AANeighborCount::GetDefaultThresholdLowHigh().GetMax()
              )
            );

            // PushBack the difference between "sl_nc" and "cb_nc" into "histogram"
            TOTAL_HISTOGRAM.PushBack( sl_nc - cb_nc);

            BCL_MessageDbg
            (
              "sl_nc is " + util::Format()( sl_nc) + " cb_nc is " + util::Format()( cb_nc) +
              " sl_nc - cb_nc is " + util::Format()( sl_nc - cb_nc)
            );
          }
          // true if "m_CalculateSLNV_CBNV" flag was set - will use neighbor vector
          else if( m_AccessibilityMeasure->GetValue() == "NeighborVector")
          {
            // create double "cb_nv" and initialize with the Norm of the neighbor vector value for the coordinate denoted
            // by "coord_itr"
            const double cb_nv
            (
              coord::NeighborVector
              (
                itr->second,
                *itr->first,
                assemble::AANeighborVector::GetDefaultThresholdLowHigh().GetMin(),
                assemble::AANeighborVector::GetDefaultThresholdLowHigh().GetMax(),
                true // normalize by neighborcount
              ).Norm()
            );

            // create double "sl_nv" and initialize with the Norm of the neighbor vector value for the spin label
            // coordinate
            const double sl_nv
            (
              coord::NeighborVector
              (
                itr->second,
                spin_label_coords,
                assemble::AANeighborVector::GetDefaultThresholdLowHigh().GetMin(),
                assemble::AANeighborVector::GetDefaultThresholdLowHigh().GetMax(),
                true // normalize by neighborcount
              ).Norm()
            );

            // add the difference between "sl_nv" and "cb_nv" to "TOTAL_HISTOGRAM"
            TOTAL_HISTOGRAM.PushBack( sl_nv - cb_nv);

            BCL_MessageDbg
            (
              "sl_nv is " + util::Format()( sl_nv) + " cb_nv is " + util::Format()( cb_nv) + " sl_nc - cb_nv is "
              + util::Format()( sl_nv - cb_nv)
            );
          }
        }
      }
    }

    //! @brief NeighborMeasureValueMeetsCriteria is for determining if a calculated neighbor measure meets a threshold
    //!        It indicates of a residue is buried enough
    //!        Important because often it is desired to only consider residues that are exposed during statistics
    //!        The neighbor measure value could be neighbor count or neighbor vector and depending on which it is
    //!        the check to see if a residue is exposed enough is inverted
    //!        i.e. (large neighbor count = buried) while (large neighbor vector = exposed)
    //! @param NEIGHBOR_MEASURE_VALUE value which will be checked for meeting the command line threshold for exposure
    //! @return bool true if "NEIGHBOR_MEASURE_VALUE" indicates the residue is exposed enough - false otherwise
    bool StatisticSpinLabel::NeighborMeasureValueMeetsCriteria( const double &NEIGHBOR_MEASURE_VALUE) const
    {
      return
      (
        m_NeighborLimitMeasure->GetValue() == "NeighborCount" &&
        NEIGHBOR_MEASURE_VALUE < m_NeighborLimit->GetNumericalValue< double>()
      )
      ||
      (
        m_NeighborLimitMeasure->GetValue() == "NeighborVector" &&
        NEIGHBOR_MEASURE_VALUE > m_NeighborLimit->GetNumericalValue< double>()
      );
    }

    //! @brief GetNeighborMeasureValue calculates a neighbor measure value for a residue based on command line measure
    //!        The neighbor measure could be neighbor count or neighbor vector; this is determined by command line
    //! @param RESIDUE_AND_NEIGHBOR_LIST the residue of interest and its neighbor list
    //! @param PRINT_ACCESSIBILITY_STATISTICS boolean indicating if each value should be printed
    //! @return double which is the neighbor measure value
    double StatisticSpinLabel::GetNeighborMeasureValue
    (
      const assemble::AANeighborList &RESIDUE_AND_NEIGHBOR_LIST,
      const bool PRINT_ACCESSIBILITY_STATISTICS
    ) const
    {
      static double accessibility_inner_product( 0);
      static double accessibility_sum( 0);
      static size_t number_accessibilities( 0);
      static double max_accessibility( 0);
      double current_count( 0);

      if( m_NeighborLimitMeasure->GetValue() == "NeighborCount")
      {
        assemble::AANeighborCount neighbor_counter;
        neighbor_counter.SetMinimalSequenceSeparation( m_NeighborLimitSeqExcl->GetNumericalValue< size_t>());

        // set "current_neighbor_count" to the neighbor count of the current amino acid denoted by "itr"
        current_count = neighbor_counter( RESIDUE_AND_NEIGHBOR_LIST);
      }
      else if( m_NeighborLimitMeasure->GetValue() == "NeighborVector")
      {
        current_count =
        assemble::AANeighborVector::NeighborVector
        (
          RESIDUE_AND_NEIGHBOR_LIST,
          assemble::AANeighborVector::GetDefaultThresholdLowHigh(),
          true // normalize by neighborcount
        ).Norm();
      }

      static double min_accessibility( current_count);
      if( !PRINT_ACCESSIBILITY_STATISTICS)
      {
        // write message giving the neighbor count/neighborvector around the amino acid currently denoted by "itr"
        BCL_MessageDbg
        (
          "the neighbor count of amino acid or neighbor vector :" +
          RESIDUE_AND_NEIGHBOR_LIST.GetCenterAminoAcid()->GetIdentification() + " is " + util::Format()( current_count)
        );
      }
      if( !PRINT_ACCESSIBILITY_STATISTICS)
      {
        accessibility_inner_product += ( current_count * current_count);
        accessibility_sum += current_count;
        ++number_accessibilities;
        if( current_count > max_accessibility)
        {
          max_accessibility = current_count;
        }
        if( current_count < min_accessibility)
        {
          min_accessibility = current_count;
        }
      }

      if( PRINT_ACCESSIBILITY_STATISTICS)
      {
        BCL_MessageCrt( "inner product " + util::Format()( accessibility_inner_product));
        BCL_MessageCrt( "accessibilities sum " + util::Format()( accessibility_sum));
        BCL_MessageCrt( "number of accessibilities " + util::Format()( number_accessibilities));
        const double accessibilites_mean( accessibility_sum / number_accessibilities);
        BCL_MessageCrt( "accessibilities mean " + util::Format()( accessibilites_mean));
        const double accessibilities_variance( math::Absolute( accessibility_inner_product / number_accessibilities - math::Sqr( accessibilites_mean)));
        BCL_MessageCrt( "accessibilities variance " + util::Format()( accessibilities_variance));
        BCL_MessageCrt( "accessibilities standard deviation " + util::Format()( math::Sqrt( accessibilities_variance)));
        BCL_MessageCrt( "minimum accessibility " + util::Format()( min_accessibility));
        BCL_MessageCrt( "maximum accessibility " + util::Format()( max_accessibility));
      }
      return current_count;
    }

    //! @brief GetAAsWithAcceptableNeighborCount gives a list of amino acids which match the neighbor count criteria
    //! @param PROTEIN_MODEL the list of all possible amino acids whose neighbor count will be checked
    //! @return returns a SiPtrVector of AABase which have acceptable neighbor countss
    const util::SiPtrVector< const biol::AABase> StatisticSpinLabel::GetAAsWithAcceptableNeighborCount
    (
      const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // create AllAANeighborList "neighbor_list" initialize with the neighbors for each amino acid
      // in "all_amino_acids"
      double distance_cutoff( assemble::AANeighborCount::GetDefaultThresholdLowHigh().GetMax());
      if( m_NeighborLimitMeasure->GetValue() == "NeighborVector")
      {
        distance_cutoff = assemble::AANeighborVector::GetDefaultThresholdLowHigh().GetMax();
      }

      const assemble::AANeighborListContainer neighbor_list
      (
        PROTEIN_MODEL.GetAminoAcids(),
        distance_cutoff,
        m_NeighborLimitSeqExcl->GetNumericalValue< size_t>(),
        true
      );

      // create SiPtrVector "exposed_amino_acids" to hold the amino acids which have an acceptable neighbor count
      util::SiPtrVector< const biol::AABase> exposed_amino_acids;

      // create double "current_neighbor_count" to hold the neighbor count/neighbor vector of each amino acid
      double current_count;

      // get amino acids which are eligable in terms of neighbor count
      for
      (
        assemble::AANeighborListContainer::const_iterator
          itr( neighbor_list.Begin()), itr_end( neighbor_list.End());
        itr != itr_end;
        ++itr
      )
      {
        current_count = GetNeighborMeasureValue( itr->second);
        // true if "current_neighbor_count" is low enough (neighbor count) or large enough (neighbor vector)
        if( NeighborMeasureValueMeetsCriteria( current_count))
        {
          BCL_MessageDbg( "adding " + itr->second.GetCenterAminoAcid()->GetIdentification() + " to exposed_amino_acids");

          // add the amino acid currently denoted by "itr" to "exposed_amino_acids"
          exposed_amino_acids.PushBack( itr->second.GetCenterAminoAcid());
        }
      }

      // return the list of amino acids with acceptable exposure
      return exposed_amino_acids;
    }

    //! @brief AddSpinLabelsAndGiveDistanceBetween calculates the coordinates of two spin labels based on the CA and
    //! CB atoms of the amino acids provided and the center of the protein and a random CA->CB->SL angle and a random
    //! CENTER->CA->CB->SL dihedral angle
    //! @param AMINO_ACID_A biol::AABackBone for which a spin label is desired to be added
    //! @param AMINO_ACID_B biol::AABackBone for which a spin label is desired to be added
    //! @param PROTEIN the protein model for which the spin labels are being added
    //! @return returns the distance between the spin labels added to "AMINO_ACID_A" and "AMINO_ACID_B"
    double StatisticSpinLabel::AddSpinLabelsAndGiveDistanceBetween
    (
      const biol::AABase &AMINO_ACID_A,
      const biol::AABase &AMINO_ACID_B,
      const assemble::ProteinModel &PROTEIN
    ) const
    {
      // create Vector3D "sl_a"
      const linal::Vector3D sl_a
      (
        GetSpinLabelCoordinates
        (
          AMINO_ACID_A,
          PROTEIN,
          m_SpinLabelLength->GetFirstParameter()->GetNumericalValue< double>(),
          m_MinCACBSLAngle->GetFirstParameter()->GetNumericalValue< double>(),
          m_MaxCACBSLAngle->GetFirstParameter()->GetNumericalValue< double>(),
          m_MinSLDihedralAngle->GetFirstParameter()->GetNumericalValue< double>(),
          m_MaxSLDihedralAngle->GetFirstParameter()->GetNumericalValue< double>()
        ).Second()
      );

      // create Vector3D "sl_b"
      const linal::Vector3D sl_b
      (
        GetSpinLabelCoordinates
        (
          AMINO_ACID_B,
          PROTEIN,
          m_SpinLabelLength->GetFirstParameter()->GetNumericalValue< double>(),
          m_MinCACBSLAngle->GetFirstParameter()->GetNumericalValue< double>(),
          m_MaxCACBSLAngle->GetFirstParameter()->GetNumericalValue< double>(),
          m_MinSLDihedralAngle->GetFirstParameter()->GetNumericalValue< double>(),
          m_MaxSLDihedralAngle->GetFirstParameter()->GetNumericalValue< double>()
        ).Second()
      );

      const double distance( linal::Distance( sl_a, sl_b));

      BCL_MessageDbg
      (
        "The distance between sl_a with coordinates : \n" + util::Format()( sl_a) + " and sl_b with coordinates \n" +
        util::Format()( sl_b) + " the distance between sl_a and sl_b is " + util::Format()( distance)
      );

      // return distance between spin labels "sl_a" and "sl_b"
      return distance;

    } // end AddSpinLabelsAndGiveDistanceBetween

    //! @brief GetSpinLabelCoordinates calculates the coordinates of s spin label based on the CA and
    //! CB atoms of the amino acid provided and the center of the protein and a random CA->CB->SL angle and a random
    //! CENTER->CA->CB->SL dihedral angle
    //! @param AMINO_ACID biol::AABackBone for which a spin label is desired to be added
    //! @param PROTEIN the protein model for which the spin labels are being added
    //! @param SL_LENGTH the length to be used for the spin label
    //! @param MIN_CA_CB_SL_ANGLE the minimum possible value for the CA->CB->SL angle
    //! @param MAX_CA_CB_SL_ANGLE the maximum possible value for the CA->CB->SL angle
    //! @param MIN_SL_DIHEDRAL_ANGLE the minimum possible value for the spin label dihedral angle
    //! @param MAX_SL_DIHEDRAL_ANGLE the minimum possible value for the spin label dihedral angle
    //! @return returns the coordinates of the spin label added to "AMINO_ACID" and the angles used to make it
    //!         The first double in the VectorND is the CA->CB->SL angle.
    //!         The second double in the VectorND is the ProtCent->CA->CB->SL dihedral angle.
    const storage::Pair< storage::VectorND< 2, double>, linal::Vector3D> StatisticSpinLabel::GetSpinLabelCoordinates
    (
      const biol::AABase &AMINO_ACID,
      const assemble::ProteinModel &PROTEIN,
      const double &SL_LENGTH,
      const double &MIN_CA_CB_SL_ANGLE,
      const double &MAX_CA_CB_SL_ANGLE,
      const double &MIN_SL_DIHEDRAL_ANGLE,
      const double &MAX_SL_DIHEDRAL_ANGLE
    ) const
    {
      // create Vector3D "amino_acid_a_cb" to hold the coordinates of the CB atom of AMINO_ACID_A
      const linal::Vector3D amino_acid_cb( AMINO_ACID.GetFirstSidechainAtom().GetCoordinates());

      BCL_Assert
      (
        AMINO_ACID.GetFirstSidechainAtom().AllCoordinatesDefined(),
        "For residue in chain " + util::Format()( AMINO_ACID.GetChainID()) + " PDB ID " +
        util::Format()( AMINO_ACID.GetPdbID()) + " all CB atom coordinates are not defined"
        + util::Format()( AMINO_ACID.GetFirstSidechainAtom().GetCoordinates())
      );

      // create Vector3D "amino_acid_a_ca" to hold the coordinates of the CA atom of AMINO_ACID_A
      const linal::Vector3D amino_acid_ca( AMINO_ACID.GetCA().GetCoordinates());

      BCL_Assert
      (
        AMINO_ACID.GetCA().AllCoordinatesDefined(),
        "For residue in chain " + util::Format()( AMINO_ACID.GetChainID()) + " PDB ID " +
        util::Format()( AMINO_ACID.GetPdbID()) + " all CA atom coordinates are not defined"
        + util::Format()( AMINO_ACID.GetCA().GetCoordinates())
      );

      // create Vector3D "protein_center" to hold the coordinates of the center of "PROTEIN"
      const linal::Vector3D protein_center( PROTEIN.GetCenter());

      // create double "ca_cb_sl_angle"
      // initialize with random angle between "minimum_ca_cb_sl_angle" and "maximum_ca_cb_sl_angle"
      double ca_cb_sl_angle( random::GetGlobalRandom().Random< double>( MIN_CA_CB_SL_ANGLE, MAX_CA_CB_SL_ANGLE));

      // create double "sl_dihedral_angle"
      // initialize with random angle between "minimum_sl_dihedral_angle" and "maximum_sl_dihedral_angle"
      double sl_dihedral_angle( random::GetGlobalRandom().Random< double>( MIN_SL_DIHEDRAL_ANGLE, MAX_SL_DIHEDRAL_ANGLE));

      // create const Vector3D "sl_coordinates"
      linal::Vector3D sl_coordinates
      (
        // initialize with the coordinates determined by "amino_acid_cb", "amino_acid_ca",
        // "sl_length", "ca_cb_sl_angle", and "sl_dihedral_angle"
        linal::CoordinatesDihedral
        (
          amino_acid_cb, amino_acid_ca, protein_center, SL_LENGTH, ca_cb_sl_angle, sl_dihedral_angle
        )
      );

      if( m_NoSLClashing->GetFlag())
      {
        const double clash_threshold( m_ClashThreshold->GetNumericalValue< double>());

        if( !SLClash( sl_coordinates, PROTEIN, clash_threshold))
        {
          // return "sl_coordinates" which is the coordinates of the spin label attached to "AMINO_ACID"
          return storage::Pair< storage::VectorND< 2, double>, linal::Vector3D>
          (
            storage::VectorND< 2, double>( ca_cb_sl_angle, sl_dihedral_angle),
            sl_coordinates
          );
        }
        const size_t max_trials( m_ClashMaxTrials->GetNumericalValue< double>());
        size_t current_trials( 0);
        bool clashing( true);
        while( clashing && current_trials < max_trials)
        {
          ca_cb_sl_angle = random::GetGlobalRandom().Random< double>( MIN_CA_CB_SL_ANGLE, MAX_CA_CB_SL_ANGLE);
          sl_dihedral_angle = random::GetGlobalRandom().Random< double>( MIN_SL_DIHEDRAL_ANGLE, MAX_SL_DIHEDRAL_ANGLE);
          sl_coordinates = linal::CoordinatesDihedral
          (
            amino_acid_cb, amino_acid_ca, protein_center, SL_LENGTH, ca_cb_sl_angle, sl_dihedral_angle
          );
          clashing = SLClash( sl_coordinates, PROTEIN, clash_threshold);
          ++current_trials;
        }

        if( clashing)
        {
          BCL_MessageTop( "failed to find nonclashing spin label");
        }
      }

      BCL_Assert
      (
        math::EqualWithinTolerance( linal::Distance( amino_acid_cb, sl_coordinates), SL_LENGTH),
        "For residue in chain " + util::Format()( AMINO_ACID.GetChainID()) + " PDB ID " +
        util::Format()( AMINO_ACID.GetPdbID()) + " distance between cb atom and sl is " +
        util::Format()( linal::Distance( amino_acid_cb, sl_coordinates)) + " but should be " + util::Format()( SL_LENGTH)
      );

      // return "sl_coordinates" which is the coordinates of the spin label attached to "AMINO_ACID"
      return storage::Pair< storage::VectorND< 2, double>, linal::Vector3D>
      (
        storage::VectorND< 2, double>( ca_cb_sl_angle, sl_dihedral_angle),
        sl_coordinates
      );

    } // end GetSpinLabelCoordinates

    //! @brief SLClash is for determining if a psuedo-spin label clashes with other atoms
    //! @param SL_COORDINATES the coordinates of the psuedo-spin label
    //! @param PROTEIN the protein in which the psuedo-spin label was placed
    //! @param CLASH_THRESHOLD the distance between the psuedo-spin label and other atoms which qualifies as a clash
    //! @return bool true if there is a clash - false if the psuedo spin label does not clash with other atoms
    bool StatisticSpinLabel::SLClash
    (
      const linal::Vector3D &SL_COORDINATES, const assemble::ProteinModel &PROTEIN, const double CLASH_THRESHOLD
    ) const
    {
      // create SiPtrVector "protein_atom_coordinates" and initialize with all the atom coordinates of "PROTEIN"
      const util::SiPtrVector< const linal::Vector3D> protein_atom_coordinates( PROTEIN.GetAtomCoordinates());

      // create bool "clash" to keep track if there is a clash between "SL_COORDINATES" and any atoms in "PROTEIN"
      // and initialize with false
      bool clash( false);

      // iterate through the atom coordinates of "protein_atom_coordinates" in order to see if any of them clash
      // with "SL_COORDINATES"
      for
      (
        util::SiPtrVector< const linal::Vector3D>::const_iterator
          atom_itr( protein_atom_coordinates.Begin()), atom_itr_end( protein_atom_coordinates.End());
        atom_itr != atom_itr_end;
        ++atom_itr
      )
      {
        // create const double "distance" and initialize with the distance between "SL_COORDINATES" and the current
        // atom coordinate denoted by "atom_itr"
        const double distance( linal::Distance( **atom_itr, SL_COORDINATES));

        // true if "distance" is less than "CLASH_THRESHOLD" - indicates there is a clash
        if( distance < CLASH_THRESHOLD)
        {
          // set "clash" to true - indicating that there is a clash
          clash = true;
        }
      }

      // return "clash"
      return clash;
    }

    //! @brief CalculateRadiusOfGyration calculate the radius of gyration of a protein model
    //! @param PROTEIN_MODEL is the protein model whose radius of gyration is desired to be calculated
    //! @return return double which is the radius of gyration
    double StatisticSpinLabel::CalculateRadiusOfGyration( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // create Set "bb_atoms" to hold the atoms desired to be included in the calculation
      storage::Set< biol::AtomType> bb_atoms;

      // loop over "ATOM_LIST_FLAG" in order to add all atoms desired to be included in the calculation to "bb_atoms"
      for
      (
        util::ShPtrVector< command::ParameterInterface>::const_iterator
          atom_itr( m_CalculateRadiusOfGyration->GetParameterList().Begin()),
          atom_itr_end( m_CalculateRadiusOfGyration->GetParameterList().End());
        atom_itr != atom_itr_end; ++atom_itr
      )
      {
        // insert the atom type into "bb_atoms"
        bb_atoms.Insert( biol::GetAtomTypes().TypeFromPDBAtomName( ( *atom_itr)->GetValue()));
      }

      BCL_MessageDbg( "The number of atoms used is : " + util::Format()( bb_atoms.GetSize()));

      // create SiPtrVector of Vector3Ds "atom_coordinates" and initialize with the atom coordinates of "PROTEIN_MODEL"
      const util::SiPtrVector< const linal::Vector3D> atom_coordinates( PROTEIN_MODEL.GetAtomCoordinates( bb_atoms));

      BCL_MessageDbg
      (
        "The number of atom coordinates is : " + util::Format()( atom_coordinates.GetSize())
      );

      // calculate the Radius of gyration using the atoms in "bb_atoms"
      return coord::RadiusOfGyration( atom_coordinates);
    }

    //! @brief function to calculate a value for SL->H - CB->H
    //! @param AA_BASE_A the residue whose amide hydrogen will be used in the current statistic
    //! @param AA_BASE_B the residue whose cb will be used in the current statistic and will have pseudo SL added
    //! @param PROTEIN_MODEL mode where the residues are located
    //! @return double which gives the value for SL->H - CB->H
    double StatisticSpinLabel::CalculateSL_NHDistance
    (
      const biol::AABase &AA_BASE_A, const biol::AABase &AA_BASE_B, const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // get the coordinates of the amide hydrogen for AA_BASE_A
      const linal::Vector3D &aa_base_a_h_coords( AA_BASE_A.GetAtom( biol::GetAtomTypes().H).GetCoordinates());

      // true if the coordinates are undefined
      if( !aa_base_a_h_coords.IsDefined())
      {
        // return undefined double
        return util::GetUndefinedDouble();
      }

      // get the distance between the amide hydrogen on aa_a and the CB atom on aa_b
      const double h_cb_distance
      (
        linal::Distance( aa_base_a_h_coords, AA_BASE_B.GetFirstSidechainAtom().GetCoordinates())
      );

      // true if the h->CB distance is undefined
      if( !util::IsDefined( h_cb_distance))
      {
        BCL_MessageDbg
        (
          "could not calculate distance between " + util::Format()( aa_base_a_h_coords) +
          " and first side chain atom of " + AA_BASE_B.GetIdentification()
        )
        // return undefined double
        return util::GetUndefinedDouble();
      }

      // get the distance between the spin label on AABASE_B and the H on AABASE_A
      const double sl_distance
      (
        linal::Distance
        (
          GetSpinLabelCoordinates
          (
            AA_BASE_B,
            PROTEIN_MODEL,
            m_SpinLabelLength->GetFirstParameter()->GetNumericalValue< double>(),
            m_MinCACBSLAngle->GetFirstParameter()->GetNumericalValue< double>(),
            m_MaxCACBSLAngle->GetFirstParameter()->GetNumericalValue< double>(),
            m_MinSLDihedralAngle->GetFirstParameter()->GetNumericalValue< double>(),
            m_MaxSLDihedralAngle->GetFirstParameter()->GetNumericalValue< double>()
          ).Second(),
          aa_base_a_h_coords
        )
      );

//      BCL_MessageStd( util::Format()(         GetSpinLabelCoordinates
//          (
//            AA_BASE_B,
//            PROTEIN_MODEL,
//            m_SpinLabelLength->GetFirstParameter()->GetNumericalValue< double>(),
//            m_MinCACBSLAngle->GetFirstParameter()->GetNumericalValue< double>(),
//            m_MaxCACBSLAngle->GetFirstParameter()->GetNumericalValue< double>(),
//            m_MinSLDihedralAngle->GetFirstParameter()->GetNumericalValue< double>(),
//            m_MaxSLDihedralAngle->GetFirstParameter()->GetNumericalValue< double>()
//          ).Second()));

      // true if the SL->H distance is undefined
      if( !util::IsDefined( sl_distance))
      {
        BCL_MessageDbg
        (
          "could not calculate distance between pseudo spin label and amide hydrogen"
        )
        // return undefined double
        return util::GetUndefinedDouble();
      }

      // calculate the difference between the SL->H distance and the CB->H distance
      const double distance_difference( sl_distance - h_cb_distance);

//      BCL_MessageStd( util::Format()( sl_distance) + "\t" + util::Format()( h_cb_distance) + "\t" + util::Format()( distance_difference));

      // return the distance difference
      return distance_difference;
    }

    //! @brief CalculateSL_CBDistance takes two AABases and returns the Dsl-Dcb value for them
    //! @param AA_BASE_A is the first AABase involved in the calculation
    //! @param AA_BASE_B is the second AABase involved in the calculation
    //! param PROTEIN_MODEL is the protein model in which "AA_BASE_A" and "AA_BASE_B" are
    //! @return double which is the Dsl-Dcb value for "AA_BASE_A" and "AA_BASE_B"
    double StatisticSpinLabel::CalculateSL_CBDistance
    (
      const biol::AABase &AA_BASE_A, const biol::AABase &AA_BASE_B, const assemble::ProteinModel &PROTEIN_MODEL
    ) const
    {
      // create const double "distance_between_spin_labels" and initialize with the distance between spin labels
      const double distance_between_spin_labels
      (
        AddSpinLabelsAndGiveDistanceBetween
        (
          AA_BASE_A,
          AA_BASE_B,
          PROTEIN_MODEL
        )
      );

      // create const double "distance_between_cb" and initialize with the distance between CB atoms
      const double distance_between_cb
      (
        linal::Distance
        (
          AA_BASE_A.GetFirstSidechainAtom().GetCoordinates(),
          AA_BASE_B.GetFirstSidechainAtom().GetCoordinates()
        )
      );

      const double dsl_dcb( distance_between_spin_labels - distance_between_cb);

      BCL_MessageDbg
      (
        "resi_a " + util::Format()( AA_BASE_A.GetPdbID()) + " chain_a " + util::Format()( AA_BASE_A.GetChainID()) +
        " resi_b " + util::Format()( AA_BASE_B.GetPdbID()) + " chain_b " + util::Format()( AA_BASE_B.GetChainID()) +
        " spin label distance is : " + util::Format()( distance_between_spin_labels) +
        " CB atoms distance " + util::Format()( distance_between_cb) + " dsl-dcb : " + util::Format()( dsl_dcb)
      );

      BCL_Assert
      (
        dsl_dcb <= 2 * m_SpinLabelLength->GetFirstParameter()->GetNumericalValue< double>() &&
        dsl_dcb >= -2 * m_SpinLabelLength->GetFirstParameter()->GetNumericalValue< double>(),
        "dsl_dcb should be between " +
        util::Format()( 2 * m_SpinLabelLength->GetFirstParameter()->GetNumericalValue< double>()) +
        " and " +
        util::Format()( -2 * m_SpinLabelLength->GetFirstParameter()->GetNumericalValue< double>()) +
        " but instead has a value of " + util::Format()( dsl_dcb)
      );

      return dsl_dcb;
    }

    //! @brief GetSSETypesFromCommandLine converts the SSE strings given over the command line into a Set of SSTypes
    //! @param FLAG is the command::FlagInterface which has the list of SSE strings
    //! @return Set of biol::SSTypes which were in the "FLAG"
    const storage::Set< biol::SSType> StatisticSpinLabel::GetSSETypesFromCommandLine( const command::FlagInterface &FLAG) const
    {
      storage::Set< biol::SSType> sse_types;

      for
      (
        util::ShPtrVector< command::ParameterInterface>::const_iterator
          param_itr( FLAG.GetParameterList().Begin()), param_itr_end( FLAG.GetParameterList().End());
        param_itr != param_itr_end;
        ++param_itr
      )
      {
        sse_types.Insert( biol::GetSSTypes().GetEnumFromName( ( *param_itr)->GetValue()));
      }
      return sse_types;
    }

    //! @brief SL_CB_SL_CBStatisticsSSESpecific does (SL-CB)-(SL-CB) statistics using specific SSE types
    //! @param PROTEIN_MODEL the protein model for which the (SL-CB)-(SL-CB) statistics will be done
    //! @param OFFSET_HISTOGRAM_VECTOR the Vector which contains the histograms for each step
    //! @param TOTAL_HISTOGRAM the histogram which all the statistics for AA_LIST will be inserted
    void StatisticSpinLabel::SL_CB_SL_CBStatisticsSSESpecific
    (
      const assemble::ProteinModel &PROTEIN_MODEL,
      storage::Vector< math::Histogram> &OFFSET_HISTOGRAM_VECTOR,
      math::Histogram &TOTAL_HISTOGRAM,
      storage::Vector< storage::VectorND< 2, storage::Vector< double> > > &DSL_VS_DCB_POINTS
    ) const
    {
      // get sse types desired to be involved in the constant spin label
      const storage::Set< biol::SSType> sse_types_constant( GetSSETypesFromCommandLine( *m_SSETypes_Constant));

      // get sse types desired to be involved in the stepping spin label
      const storage::Set< biol::SSType> sse_types_step( GetSSETypesFromCommandLine( *m_SSETypes_Step));

      // create SiPtrVector "all_amino_acids" and initialize with the amino acids in the SSEs of "protein_model"
      const util::SiPtrVector< const biol::AABase> all_amino_acids( PROTEIN_MODEL.GetAminoAcids());

      static const double s_max_distance( 1000.0);

      // create AllAANeighborList "neighbor_list" initialize with the neighbors for each amino acid
      // in "all_amino_acids"
      const assemble::AANeighborListContainer neighbor_list
      (
        all_amino_acids, s_max_distance, 0, true
      );

      // true if all SSE types are desired to be used so can just iterate straight over all the residues of the protein
      if
      (
        sse_types_constant.GetSize() == biol::GetSSTypes().COIL.GetIndex() + 1 &&
        sse_types_step.GetSize() == biol::GetSSTypes().COIL.GetIndex() + 1
      )
      {
        CalculateSL_CB_SL_CBStatistics
        (
          all_amino_acids,         all_amino_acids, neighbor_list,
          OFFSET_HISTOGRAM_VECTOR, TOTAL_HISTOGRAM, PROTEIN_MODEL, DSL_VS_DCB_POINTS
        );
      }

      // else true if a subset of SSEs for the constant and the stepping residues are specified as a desired to be used
      // or also true if only a subset of SSEs for the stepping residues are specified as a desired to be used
      // Need to iterate over the AASequence of each SSE rather than straight through the whole protein's AASequence.
      // It is okay to break up the AASequence of the protein into SSE AASequences for the constant residue (even
      // though it might not be necessary, i.e. in the second part of this if statement)
      // since all that is cared about for the constant residue is that it somehow
      // goes on every allowable residue. For the stepping residue, we don't want to break up the protein sequence
      // into the smaller SSE AASequences if we don't need to because this limits the walking distance that can be
      // achieved.
      else if
      (
        (
          sse_types_constant.GetSize() < biol::GetSSTypes().COIL.GetIndex() + 1 &&
          sse_types_step.GetSize() < biol::GetSSTypes().COIL.GetIndex() + 1
        )
        ||
        (
          sse_types_constant.GetSize() == biol::GetSSTypes().COIL.GetIndex() + 1 &&
          sse_types_step.GetSize() < biol::GetSSTypes().COIL.GetIndex() + 1
        )
      )
      {
        // create SiPtr Vector "sses_constant" and initialize with the sses which are desired to be
        // used with the constant residue
        util::SiPtrVector< const assemble::SSE> sses_constant( PROTEIN_MODEL.GetSSEs( sse_types_constant));

        // create SiPtr Vector "sses_step" and initialize with the sses which are desired to be
        // used with the stepping residue
        util::SiPtrVector< const assemble::SSE> sses_step( PROTEIN_MODEL.GetSSEs( sse_types_step));

        // iterate over all of the sses in "sses_constant"
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator
            sse_list_itr_a( sses_constant.Begin()), sse_list_itr_end( sses_constant.End());
          sse_list_itr_a != sse_list_itr_end;
          ++sse_list_itr_a
        )
        {
          // get the amino acids in the sse currently denoted by "sse_list_itr_a"
          const util::ShPtrVector< biol::AABase> sse_aa_list_a( ( *sse_list_itr_a)->GetData());

          // iterate over all of the sses in "sses_step"
          for
          (
            util::SiPtrVector< const assemble::SSE>::const_iterator
              sse_list_itr_b( sses_step.Begin()), sse_list_itr_b_end( sses_step.End());
            sse_list_itr_b != sse_list_itr_b_end;
            ++sse_list_itr_b
          )
          {
            // get the amino acids in the sse currently denoted by "sse_list_itr_b"
            const util::ShPtrVector< biol::AABase> sse_aa_list_b( ( *sse_list_itr_b)->GetData());

            // calculate the statistics for the aas in "sse_aa_list_a" and "sse_aa_list_b"
            CalculateSL_CB_SL_CBStatistics
            (
              sse_aa_list_a,           sse_aa_list_b,   neighbor_list,
              OFFSET_HISTOGRAM_VECTOR, TOTAL_HISTOGRAM, PROTEIN_MODEL, DSL_VS_DCB_POINTS
            );
          }
        }
      }

      // else true if only a subset of SSEs for the constant residues are specified as a desired to be used
      else if
      (
        sse_types_constant.GetSize() < biol::GetSSTypes().COIL.GetIndex() + 1 &&
        sse_types_step.GetSize() == biol::GetSSTypes().COIL.GetIndex() + 1
      )
      {
        // create SiPtr Vector "sses_constant" and initialize with the sses which are desired to be
        // used with the constant residue
        util::SiPtrVector< const assemble::SSE> sses_constant( PROTEIN_MODEL.GetSSEs( sse_types_constant));

        // create SiPtrVector "all_amino_acids" and initialize with the amino acids in the SSEs of "protein_model"
        const util::SiPtrVector< const biol::AABase> all_amino_acids( PROTEIN_MODEL.GetAminoAcids());

        // iterate over all of the sses in "sses_constant"
        for
        (
          util::SiPtrVector< const assemble::SSE>::const_iterator
            sse_list_itr_a( sses_constant.Begin()), sse_list_itr_end( sses_constant.End());
          sse_list_itr_a != sse_list_itr_end;
          ++sse_list_itr_a
        )
        {
          // get the amino acids in the sse currently denoted by "sse_list_itr_a"
          const util::ShPtrVector< biol::AABase> sse_aa_list_a( ( *sse_list_itr_a)->GetData());

          // calculate the statistics for the aas in "sse_aa_list_a" and "sse_aa_list_b"
          CalculateSL_CB_SL_CBStatistics
          (
            sse_aa_list_a,           all_amino_acids,   neighbor_list,
            OFFSET_HISTOGRAM_VECTOR, TOTAL_HISTOGRAM, PROTEIN_MODEL, DSL_VS_DCB_POINTS
          );
        }
      }
      else
      {
        BCL_Exit( "In SL_CB_SL_CBStatisticsSSESpecific, could not determine which residues to iterate over.", -1);
      }
    }

    //! @brief CalculateSL_CB_SL_CBStatistics does statistics for difference between two SL-CB measurements
    //! (Dsl-Dcb)-(Dsl-Dcb)
    //! @param AA_LIST_A TODO: document
    //! @param AA_LIST_B TODO: document
    //! @param NEIGHBOR_LIST TODO: document
    //! @param OFFSET_HISTOGRAM_VECTOR the Vector which contains the histograms for each step
    //! @param TOTAL_HISTOGRAM the histogram which all the statistics for AA_LIST will be inserted
    //! @param PROTEIN_MODEL is the protein model which "AA_LIST" corresponds to
    //! @param DSL_VS_DCB_POINTS TODO: document
    template< typename t_DataType_A, typename t_DataType_B>
    void StatisticSpinLabel::CalculateSL_CB_SL_CBStatistics
    (
      const storage::Vector< t_DataType_A> &AA_LIST_A,
      const storage::Vector< t_DataType_B> &AA_LIST_B,
      const assemble::AANeighborListContainer &NEIGHBOR_LIST,
      storage::Vector< math::Histogram> &OFFSET_HISTOGRAM_VECTOR,
      math::Histogram &TOTAL_HISTOGRAM,
      const assemble::ProteinModel &PROTEIN_MODEL,
      storage::Vector< storage::VectorND< 2, storage::Vector< double> > > &DSL_VS_DCB_POINTS
    ) const
    {
      // loop over amino acids in "AA_LIST" to get the (Dsl-Dcb)-(Dsl-Dcb) distances
      for
      (
        typename storage::Vector< t_DataType_A>::const_iterator
          aa_itr_a( AA_LIST_A.Begin()), aa_itr_a_end( AA_LIST_A.End());
        aa_itr_a != aa_itr_a_end;
        ++aa_itr_a
      )
      {
        // true if any of the coordinates are not defined for the CB of the anchor aa currently denoted by "aa_itr_a"
        if( !( *aa_itr_a)->GetFirstSidechainAtom().AllCoordinatesDefined() || !( *aa_itr_a)->GetCA().AllCoordinatesDefined())
        {
          BCL_MessageDbg
          (
            "CB or CA atom coordinates not defined for residue : \n" + util::Format()( **aa_itr_a)
          );
          // coordinates are not defined so go to next amino acid for anchor residue
          continue;
        }
        // this if statement needs to come after the previous if statement because if the coordinates aren't defined
        // then they residue won't be anywhere in "NEIGHBOR_LIST" and this will trigger the following assert
        if( m_LimitNeighbors->GetFlag())
        {
          BCL_Assert
          (
            NEIGHBOR_LIST.Find( **aa_itr_a) != NEIGHBOR_LIST.End(),
            "aa_itr_a not found in NEIGHBOR_LIST " + util::Format()( **aa_itr_a)
          );
          if( !NeighborMeasureValueMeetsCriteria( GetNeighborMeasureValue( NEIGHBOR_LIST.Find( **aa_itr_a)->second)))
          {
            continue;
          }
        }

        // loop over amino acids in "AA_LIST" to get the (Dsl-Dcb)-(Dsl-Dcb) statistics
        for
        (
          typename storage::Vector< t_DataType_B>::const_iterator
            aa_itr_b( AA_LIST_B.Begin()), aa_itr_b_end( AA_LIST_B.End());
          aa_itr_b != aa_itr_b_end;
          ++aa_itr_b
        )
        {

          // true if any coordinates are not defined for the CB of the reference aa currently denoted by "aa_itr_b"
          // or true if "aa_itr_b" is pointing to the same aa as "aa_itr_a"
          if
          (
            !( *aa_itr_b)->GetFirstSidechainAtom().AllCoordinatesDefined() || !( *aa_itr_b)->GetCA().AllCoordinatesDefined() ||
            *aa_itr_b == *aa_itr_a
          )
          {
            BCL_MessageDbg
            (
              "CB or CA atom coordinates not defined for residue : \n" + util::Format()( **aa_itr_b)
            );
            // coordinates are not defined or "aa_itr_b" and "aa_itr_a" point to the same aa so go to next amino acid
            // for reference residue
            continue;
          }

          // this if statement needs to come after the previous if statement because if the coordinates aren't defined
          // then they residue won't be anywhere in "NEIGHBOR_LIST" and this will trigger the following assert
          if( m_LimitNeighbors->GetFlag())
          {
            BCL_Assert
            (
              NEIGHBOR_LIST.Find( **aa_itr_b) != NEIGHBOR_LIST.End(),
              "aa_itr_b not found in NEIGHBOR_LIST " + util::Format()( **aa_itr_b)
            )
            if
            ( !NeighborMeasureValueMeetsCriteria( GetNeighborMeasureValue( NEIGHBOR_LIST.Find( **aa_itr_b)->second)))
            {
              continue;
            }
          }

          const command::ParameterInterface &spin_label_length_param( *m_SpinLabelLength->GetFirstParameter());
          const command::ParameterInterface &min_cacb_sl_angle_param( *m_MinCACBSLAngle->GetFirstParameter());
          const command::ParameterInterface &max_cacb_sl_angle_param( *m_MaxCACBSLAngle->GetFirstParameter());
          const command::ParameterInterface &min_sl_dihedral_angle_param( *m_MinSLDihedralAngle->GetFirstParameter());
          const command::ParameterInterface &max_sl_dihedral_angle_param( *m_MaxSLDihedralAngle->GetFirstParameter());

          // create "const_sl_coords" to hold coords of the spin label which will not move i.e. on the anchor residue
          const storage::Pair< storage::VectorND< 2, double>, linal::Vector3D> const_sl_coords
          (
            GetSpinLabelCoordinates
            (
              **aa_itr_a,
              PROTEIN_MODEL,
              spin_label_length_param.GetNumericalValue< double>(),
              min_cacb_sl_angle_param.GetNumericalValue< double>(),
              max_cacb_sl_angle_param.GetNumericalValue< double>(),
              min_sl_dihedral_angle_param.GetNumericalValue< double>(),
              max_sl_dihedral_angle_param.GetNumericalValue< double>()
            )
          );

          // create "const_cb_coords" to hold coords of the CB which will not move i.e. on the anchor residue
          const linal::Vector3D const_cb_coords( ( *aa_itr_a)->GetFirstSidechainAtom().GetCoordinates());

          // create "initial_sl_coords" to hold coords of the spin label which will move i.e. on the walking residue
          const storage::Pair< storage::VectorND< 2, double>, linal::Vector3D> initial_sl_coords
          (
            GetSpinLabelCoordinates
            (
              **aa_itr_b,
              PROTEIN_MODEL,
              spin_label_length_param.GetNumericalValue< double>(),
              min_cacb_sl_angle_param.GetNumericalValue< double>(),
              max_cacb_sl_angle_param.GetNumericalValue< double>(),
              min_sl_dihedral_angle_param.GetNumericalValue< double>(),
              max_sl_dihedral_angle_param.GetNumericalValue< double>()
            )
          );

          // create "initial_cb_coords" to hold coords of the CB which will move i.e. on the walking residue
          const linal::Vector3D initial_cb_coords( ( *aa_itr_b)->GetFirstSidechainAtom().GetCoordinates());

          // create double "dsl_initial" initialize with the distance between "const_sl_coords" and "initial_sl_coords"
          const double dsl_initial( linal::Distance( const_sl_coords.Second(), initial_sl_coords.Second()));

          // create double "dcb_initial" initialize with the distance between "const_cb_coords" and "initial_cb_coords"
          const double dcb_initial( linal::Distance( const_cb_coords, initial_cb_coords));

          // create const double "sl_cb_a" and initialize with the Dsl-Dcb value
          const double sl_cb_a( dsl_initial - dcb_initial);

          const int aa_itr_b_position( aa_itr_b - AA_LIST_B.Begin());

          // create int "desired_start_position" and initialize with the index to start moving spin label position
          const int desired_start_position
          (
            aa_itr_b_position +
            m_StartOffset->GetNumericalValue< int>()
          );

          // create int "start_position" and initialize with "desired_start_position"
          int start_position( desired_start_position);

          // true if "desired_start_position" is less than zero
          if( desired_start_position < 0)
          {
            // can't access negative position in "AA_LIST" so set "start_position" to zero i.e. shift to zero
            start_position = 0;
          }

          // move the spin label a given number of positions and calculate (SL-CB)-(SL-CB) for each position
          // hold "aa_itr_a" constant
          for
          (
            int position( start_position),
                // take into account to reduce the number of steps if "desired_start_position" had to be shifted
                end_position
                (
                    start_position
                  + m_NumberOFSteps->GetNumericalValue< int>()
                  - ( start_position - desired_start_position)
                );
            position < end_position && position < int( AA_LIST_B.GetSize()); //< make sure don't go past end of AAlist
            ++position
          )
          {
            // create int "aa_itr_a_position" and initialize with the position of "aa_itr" in "AA_LIST_A"
            const int aa_itr_a_position( aa_itr_a - AA_LIST_A.Begin());

            // create size_t "histogram_number" initialize with the histogram the current calculation will be put into
            const size_t histogram_number( position - start_position + ( start_position - desired_start_position));

            BCL_MessageDbg
            (
              "Starting position : "
              + util::Format()( start_position) + " Current position : " + util::Format()( position)
              + " desired start position : " + util::Format()( desired_start_position)
              + " reference position : " + util::Format()( aa_itr_b_position)
              + " constant position : " + util::Format()( aa_itr_a_position)
              + " histogram number : " + util::Format()( histogram_number)
            );
            // true if "position" denotes the same amino acid as "aa_itr_b"
            // or true if the coordinates CB or CA for the walking aa are not defined
            if
            (
              position == aa_itr_b_position || !AA_LIST_B( position)->GetFirstSidechainAtom().AllCoordinatesDefined() ||
              !AA_LIST_B( position)->GetCA().AllCoordinatesDefined()
            )
            {
              BCL_MessageDbg
              (
                "skipping position : " + util::Format()( position)
                + " should be same as reference position : " + util::Format()( aa_itr_b_position)
              );

              // skip this calculation
              continue;
            }
            // this if statement needs to come after the previous if statement because if the coordinates aren't defined
            // then they residue won't be anywhere in "NEIGHBOR_LIST" and this will trigger the following assert
            if( m_LimitNeighbors->GetFlag())
            {
              BCL_Assert
              (
                NEIGHBOR_LIST.Find( *AA_LIST_B( position)) != NEIGHBOR_LIST.End(),
                "position not found in NEIGHBOR_LIST " + util::Format()( *AA_LIST_B( position))
              )
              if
              (
                !NeighborMeasureValueMeetsCriteria( GetNeighborMeasureValue( NEIGHBOR_LIST.Find( *AA_LIST_B( position))->second))
              )
              {
                continue;
              }
            }

            // make sure that the walking spin label does not walk over to another chain
            if( !( AA_LIST_B( position)->GetChainID() == ( *aa_itr_b)->GetChainID()))
            {
              continue;
            }

            // create double "walking_ca_cb_sl_angle_min"
            // initialize with the CA->CB->SL angle for the initial spin label minus the command line given drift
            const double walking_ca_cb_sl_angle_min
            (
              initial_sl_coords.First().First() -
              m_CaCbSlAngleDrift->GetNumericalValue< int>()
            );

            // create double "walking_ca_cb_sl_angle_max"
            // initialize with the CA->CB->SL angle for the initial spin label plus the command line given drift
            const double walking_ca_cb_sl_angle_max
            (
              initial_sl_coords.First().First() +
              m_CaCbSlAngleDrift->GetNumericalValue< int>()
            );

            // create double "walking_sl_dihedral_angle_min" initialize with the proteincenter->CA->CB->SL dihedral
            //  angle for the initial spin label minus the command line given drift
            const double walking_sl_dihedral_angle_min
            (
              initial_sl_coords.First().Second() -
              m_SlDihedralAngleDrift->GetNumericalValue< int>()
            );

            // create double "walking_sl_dihedral_angle_max" initialize with the proteincenter->CA->CB->SL dihedral
            //  angle for the initial spin label plus the command line given drift
            const double walking_sl_dihedral_angle_max
            (
              initial_sl_coords.First().Second() +
              m_SlDihedralAngleDrift->GetNumericalValue< int>()
            );

            // create "walking_sl_coords" to hold coords of the spin label which is moving i.e. the walking residue
            const storage::Pair< storage::VectorND< 2, double>, linal::Vector3D> walking_sl_coords
            (
              GetSpinLabelCoordinates
              (
                *AA_LIST_B( position), //< amino acid to get SL for
                PROTEIN_MODEL,         //< protein model in which the amino acid is
                spin_label_length_param.GetNumericalValue< double>(), //< sl length
                walking_ca_cb_sl_angle_min,
                walking_ca_cb_sl_angle_max,
                walking_sl_dihedral_angle_min,
                walking_sl_dihedral_angle_max
              )
            );

            // create "walking_cb_coords" to hold coords of the CB which will move i.e. on the walking residue
            const linal::Vector3D walking_cb_coords( AA_LIST_B( position)->GetFirstSidechainAtom().GetCoordinates());

            // create double "dsl_walking" initialize with distance between "const_sl_coords" and "walking_sl_coords"
            const double dsl_walking( linal::Distance( const_sl_coords.Second(), walking_sl_coords.Second()));

            BCL_Assert
            (
              util::IsDefined( dsl_walking), "dsl_walking is not defined when created from coordinates " +
              util::Format()( const_sl_coords.Second()) + "and" + util::Format()( walking_sl_coords.Second())
            );

            // create double "dcb_walking" initialize with distance between "const_cb_coords" and "walking_cb_coords"
            const double dcb_walking( linal::Distance( const_cb_coords, walking_cb_coords));

            // create double "sl_cb_b" and initialize with the Dsl-Dcb value between "aa_itr_a" and the AABase at
            // "position" in "AA_LIST_B
            const double sl_cb_b( dsl_walking - dcb_walking);

            BCL_MessageDbg
            (
              "for position " + util::Format()( histogram_number) +
              " : Dsl_initial-Dsl_walking is " + util::Format()( dsl_initial - dsl_walking) +
              " : Dcb_initial-Dcb_walking is " + util::Format()( dcb_initial - dcb_walking)
            );

            BCL_MessageDbg
            (
              "for position " + util::Format()( histogram_number) +
              " : Dsl-Dcb_initial is " + util::Format()( sl_cb_a) +
              " : Dsl-Dcb_walking is " + util::Format()( sl_cb_b)
            );

            // make sure "histogram_number" is a valid position in "OFFSET_HISTOGRAM_VECTOR"
            BCL_Assert
            (
              histogram_number < OFFSET_HISTOGRAM_VECTOR.GetSize(),
              "histogram number not less than size of histogram vector : size is "
              + util::Format()( OFFSET_HISTOGRAM_VECTOR.GetSize()) + "\n Starting position : "
              + util::Format()( start_position) + " Current position : " + util::Format()( position)
              + " desired start position : " + util::Format()( desired_start_position)
              + " calculated histogram number : " + util::Format()( histogram_number)
              + " reference position : " + util::Format()( aa_itr_b_position)
            );
            BCL_Assert
            (
              util::IsDefined( dsl_initial - dsl_walking),
              "dsl_initial - dsl_walking is not defined. dsl_initial is " + util::Format()( dsl_initial) +
              " and dsl_walking is " + util::Format()( dsl_walking)
            );

            if( m_RecordDslVsDCB->GetValue() == "true")
            {
              DSL_VS_DCB_POINTS( histogram_number).First().PushBack( dsl_initial - dsl_walking);
              BCL_Assert( util::IsDefined( dcb_initial - dcb_walking), "util::IsDefined( dcb_initial - dcb_walking)");
              DSL_VS_DCB_POINTS( histogram_number).Second().PushBack( dcb_initial - dcb_walking);
            }
            // add the difference between "sl_cb_a" and "sl_cb_b" to the appropriate specific histogram
            OFFSET_HISTOGRAM_VECTOR( histogram_number).PushBack( sl_cb_a - sl_cb_b);

            // add the difference between "sl_cb_a" and "sl_cb_b" to "TOTAL_HISTOGRAM"
            TOTAL_HISTOGRAM.PushBack( sl_cb_a - sl_cb_b);
          }
        }
      }
    }

    //! @brief CreateHistogram2DFromValuePairFile is for creating a histogram 2d and gnuplot script from a file
    void StatisticSpinLabel::CreateHistogram2DFromValuePairFile() const
    {
      // read in the data from "m_Histogram2DValuePairFilename"
      io::IFStream read;
      io::File::MustOpenIFStream( read, m_Histogram2DValuePairFilename->GetValue());
      const storage::Vector< storage::Vector< std::string> > value_pairs( util::SplittedStringLineListFromIStream( read));
      io::File::CloseClearFStream( read);

      // create const storage::VectorND< 2, double> "min_x_y" and initialize with "m_MinX" and "m_MinY"
      const storage::VectorND< 2, double> min_x_y
      (
        m_Histogram2DMinX->GetNumericalValue< double>(), m_Histogram2DMinY->GetNumericalValue< double>()
      );

      // create const storage::VectorND< 2, double> "binsize_x_y" and initialize with "m_BinSizeX" and "m_BinSizeY"
      const storage::VectorND< 2, double> binsize_x_y
      (
        m_Histogram2DBinSizeX->GetNumericalValue< double>(), m_Histogram2DBinSizeY->GetNumericalValue< double>()
      );

      // create const storage::VectorND< 2, double> "num_bins_x_y" and initialize with "m_NumBinsX" and "m_NumBinsY"
      const storage::VectorND< 2, size_t> num_bins_x_y
      (
        m_Histogram2DNumBinsX->GetNumericalValue< size_t>(), m_Histogram2DNumBinsY->GetNumericalValue< size_t>()
      );

      // create histogram2d "histogram" and initialize with "min_x_y", "binsize_x_y", and "num_bins_x_y"
      math::Histogram2D histogram( min_x_y, binsize_x_y, num_bins_x_y);

      // iterate through "value_pairs" in order to fill "histogram"
      for
      (
        storage::Vector< storage::Vector< std::string> >::const_iterator
          itr( value_pairs.Begin()), itr_end( value_pairs.End());
        itr != itr_end;
        ++itr
      )
      {
        // get const reference "current_values" to vector currently denoted by "itr"
        const storage::Vector< std::string> &current_values( *itr);

        BCL_MessageDbg( "current_resi is :" + util::Format()( current_values( 0)));

        // convert the strings in "current_values" to doubles and create a VectorND and pushback it into "histogram"
        histogram.PushBack
        (
          storage::VectorND< 2, double>
          (
            util::ConvertStringToNumericalValue< double>
            (
              current_values( m_Histogram2DValueColumnX->GetNumericalValue< size_t>())
            ),
            util::ConvertStringToNumericalValue< double>
            (
              current_values( m_Histogram2DValueColumnY->GetNumericalValue< size_t>())
            )
          )
        );
      }

      if( util::ConvertStringToBoolean( m_Histogram2DNormalize->GetValue()))
      {
        // normalize "histogram"
        histogram.Normalize();
      }

      // write "histogram" to file designated by "m_Histogram2DOutputFilename"
      io::OFStream write;
      io::File::MustOpenOFStream( write, m_Histogram2DOutputFilename->GetValue());
      write << histogram;
      io::File::CloseClearFStream( write);

      // write gnuplot for "histogram to file designated by "m_Histogram2DGnuplotOutputFilename"
      io::File::MustOpenOFStream( write, m_Histogram2DGnuplotOutputFilename->GetValue());
      math::GnuplotHeatmap heatmap;
      heatmap.SetFromHistogram( histogram, true, true);
      heatmap.SetTitleAndLabel( m_Histogram2DGnuplotTitle->GetValue(), "", "", "");
      heatmap.WriteScript( write);
      io::File::CloseClearFStream( write);
    }

    //! @brief CreateHistogramFromValueFile is for creating a histogram and gnuplot script from a file of values
    void StatisticSpinLabel::CreateHistogramFromValueFile() const
    {
      // read in the data from "m_HistogramValueFilename"
      io::IFStream read;
      BCL_MessageStd( "opening " + m_HistogramValueFilename->GetValue());
      io::File::MustOpenIFStream( read, m_HistogramValueFilename->GetValue());
      const storage::Vector< storage::Vector< std::string> > values( util::SplittedStringLineListFromIStream( read));
      io::File::CloseClearFStream( read);

      // create const double "min_value" and initialize with "m_HistogramMin"
      const double min( m_HistogramMin->GetNumericalValue< double>());

      // create const double "bin_size" and initialize with "m_HistogramBinSize"
      const double bin_size( m_HistogramBinSize->GetNumericalValue< double>());

      // create const double "num_bins" and initialize with "m_HistogramNumBins"
      const double num_bins( m_HistogramNumBins->GetNumericalValue< size_t>());

      // create histogram2d "histogram" and initialize with "min", "bin_size", and "num_bins"
      math::Histogram histogram( min, bin_size, num_bins);

      // create const size_t "value_column" and initialize with "m_HistogramValueColumn"
      const size_t value_column( m_HistogramValueColumn->GetNumericalValue< size_t>());

      // create const size_t "weight_column" and initialize with "m_HistogramWeightColumn"
      const size_t weight_column( m_HistogramWeightColumn->GetNumericalValue< size_t>());

      // iterate through "values" in order to fill "histogram"
      for
      (
        storage::Vector< storage::Vector< std::string> >::const_iterator
          itr( values.Begin()), itr_end( values.End());
        itr != itr_end;
        ++itr
      )
      {
        // get const reference "current_values" to vector currently denoted by "itr"
        const storage::Vector< std::string> &current_values( *itr);

        BCL_Assert
        (
          current_values.GetSize() > value_column, "string vector has only " + util::Format()( current_values.GetSize())
          + " strings but value column is " + util::Format()( value_column) + "\nthe current string vector is "
          + util::Format()( current_values)
        );

        // create double "current_weight" and initialize with a value of 1. this will be changed below if weights are
        // actually needed as determined by if the "weight_column" is defined
        double current_weight( 1);

        // true if "weight_column" is defined
        if( util::IsDefined( weight_column))
        {
          // make sure the weight_column is within the bounds of the number of entries on the current line
          BCL_Assert
          (
            current_values.GetSize() > weight_column, "string vector has only " + util::Format()( current_values.GetSize())
            + " strings but weight_column column is " + util::Format()( weight_column) + "\nthe current string vector is "
            + util::Format()( current_values)
          );

          // set "weight_column" to the value that is provided by the current line "current_values"
          current_weight = util::ConvertStringToNumericalValue< double>( current_values( weight_column));
        }

        // convert the string in index "value_column" of "current_values" to double and do the same for the weight
        const double current_value( util::ConvertStringToNumericalValue< double>( current_values( value_column)));

        BCL_MessageDbg
        (
          "current value is " + util::Format()( current_value)
          + " current weight is " + util::Format()( current_weight)
        );

        // convert the string in index "value_column" of "current_values" to double and pushback it into "histogram"
        histogram.PushBack( current_value, current_weight);
      }

      if( util::ConvertStringToBoolean( m_HistogramNormalize->GetValue()))
      {
        // normalize "histogram"
        histogram.Normalize();
      }

      // write "histogram" to file designated by "m_Histogram2DOutputFilename"
      io::OFStream write;
      io::File::MustOpenOFStream( write, m_HistogramOutputFilename->GetValue());
      write << histogram;
      io::File::CloseClearFStream( write);

      // write gnuplot for "histogram to file designated by "m_HistogramGnuplotOutputFilename"
      io::File::MustOpenOFStream( write, m_HistogramGnuplotOutputFilename->GetValue());
      math::GnuplotHeatmap heatmap;
      heatmap.SetFromHistogram( histogram, false, true);
      heatmap.SetTitleAndLabel( m_HistogramGnuplotOutputFilename->GetValue(), "", "", "");
      heatmap.SetFilename( io::File::RemoveLastExtension( m_HistogramGnuplotOutputFilename->GetValue()));
      heatmap.WriteScript( write);
      io::File::CloseClearFStream( write);
    }

    //! @brief CalculateEPRDistanceAgreement for calculating the agreement of a protein model with distance restraints
    void StatisticSpinLabel::CalculateEPRDistanceAgreement() const
    {
      // create std::string "restraint_filename" and initialize with the name of the file containing the restraints
      std::string restraint_filename( m_DistanceRestraintFile->GetValue());

      // create io::IFStream "read"
      io::IFStream read;

      // open "read" with "restraint_filename"
      io::File::MustOpenIFStream( read, restraint_filename);

      // create ShPtr to ShPtrVector of restraints "restraints" and initialize with restraints from "read"
      util::ShPtrVector< restraint::AtomDistance> restraints
      (
        restraint::HandlerAtomDistanceAssigned().ReadRestraints( read)
      );

      // close and clear "read"
      io::File::CloseClearFStream( read);

      // table to hold the results
      storage::TableHeader header;

      // fill table header
      // iterate through the restraints
      for
      (
        util::ShPtrVector< restraint::AtomDistance>::const_iterator
          restraint_itr( restraints.Begin()), restraint_itr_end( restraints.End());
        restraint_itr != restraint_itr_end;
        ++restraint_itr
      )
      {
        header.PushBack( ( *restraint_itr)->GetIdentification());
      }

      header.PushBack( "sum");

      storage::Table< double> results_table( header);

      score::RestraintDistanceSpinLabel score;

      // read all the pdb list filenames from the provided file and set pdb_lists to the list of pdb list filenames
      const std::string &list_of_pdb_lists_filename( m_EPRDistanceAgreementPDBList->GetValue());
      io::File::MustOpenIFStream( read, list_of_pdb_lists_filename);
      storage::Vector< std::string> pdb_list_filenames( util::StringLineListFromIStream( read));
      io::File::CloseClearFStream( read);

      // one dataset for each restraint plus the scoresum
      storage::Vector< math::RunningAverageSD< double> > restraint_statistics( restraints.GetSize() + 1);

      // iterate through the pdb list filenames
      for
      (
        storage::Vector< std::string>::const_iterator
          pdb_list_filename_itr( pdb_list_filenames.Begin()), pdb_list_filename_itr_end( pdb_list_filenames.End());
        pdb_list_filename_itr != pdb_list_filename_itr_end;
        ++pdb_list_filename_itr
      )
      {
        // get current protein
        const assemble::ProteinModel model( GetProteinModel( *pdb_list_filename_itr));
        storage::Vector< double> current_scores;
        double scoresum( 0);
        // iterate through the restraints
        auto statistics_itr( restraint_statistics.Begin()), statistics_itr_end( restraint_statistics.End());
        for
        (
          util::ShPtrVector< restraint::AtomDistance>::const_iterator
            restraint_itr( restraints.Begin()), restraint_itr_end( restraints.End());
          restraint_itr != restraint_itr_end && statistics_itr != statistics_itr_end;
          ++restraint_itr, ++statistics_itr
        )
        {
          const restraint::AtomDistanceAssignment assignment( ( *restraint_itr)->GenerateAssignment( model));
          double current_score( score( assignment));
          if( util::ConvertStringToBoolean( m_CalculateDistanceDifferences->GetValue()))
          {
            current_score = assignment.CalculateAtomDistance() - assignment.GetDistance();
          }
          current_scores.PushBack( current_score);
          scoresum += current_score;
          *statistics_itr += current_score;
        }
        current_scores.PushBack( scoresum);
        results_table.InsertRow( *pdb_list_filename_itr, current_scores);
        restraint_statistics.LastElement() += scoresum;
      }

      // get the average agreement of each restraint and put into table
      {
        storage::Vector< double> statistic;
        for
        (
          auto
            statistics_itr( restraint_statistics.Begin()), statistics_itr_end( restraint_statistics.End());
          statistics_itr != statistics_itr_end;
          ++statistics_itr
         )
        {
          statistic.PushBack( statistics_itr->GetAverage());
        }
        results_table.InsertRow( "average", statistic);
      }
      // get the standard deviation of agreement of each restraint and put into table
      {
        storage::Vector< double> statistic;
        for
        (
          auto
            statistics_itr( restraint_statistics.Begin()), statistics_itr_end( restraint_statistics.End());
          statistics_itr != statistics_itr_end;
          ++statistics_itr
         )
        {
          statistic.PushBack( statistics_itr->GetStandardDeviation());
        }
        results_table.InsertRow( "stddev", statistic);
      }

      // create std::string "" and initialize with the name of the score file given over the command line
      std::string score_output_filename( m_DistanceRestraintScoreOutputFile->GetValue());

      // create IFStream write
      io::OFStream write;

      // open "write" with "score_output_filename"
      io::File::MustOpenOFStream( write, score_output_filename);

      results_table.WriteFormatted( write);

      // close and clear "write"
      io::File::CloseClearFStream( write);
    }

    //! @brief GetProteinModel is for creating a protein model given a pdb filename
    //! @param PDB_FILENAME is the name of the pdb filename from which a protein model will be created
    //! @return assemble::ProteinModel created from "PDB_FILENAME"
    assemble::ProteinModel StatisticSpinLabel::GetProteinModel( const std::string &PDB_FILENAME) const
    {
      // output of current pdb filename
      BCL_MessageStd( "processing pdb: " + PDB_FILENAME);

      // create static PDBFactory "pdb_factory" and initialize it to work with AABACKBONE type amino acids
      static const pdb::Factory pdb_factory;

      //// create io::IFStream "read"
      io::IFStream read;

      // open "read" and bind it to "pdb_filename"
      io::File::MustOpenIFStream( read, PDB_FILENAME);

      // create PDBReader "pdb" from "read"
      pdb::Handler pdb( read);

      // close and clear "read"
      io::File::CloseClearFStream( read);

      // message telling which pdb is being processed
      BCL_MessageStd( "processing pdb: " + PDB_FILENAME + " is complete");

      // create ProteinModel "protein_model" and initialize with "pdb"
      const assemble::ProteinModel protein_model( pdb_factory.ProteinModelFromPDB( pdb));

      // return "protein_model"
      return protein_model;
    }

    // default constructor
    StatisticSpinLabel::StatisticSpinLabel() :
      m_LimitNeighbors
      (
        new command::FlagStatic
        (
          "limit_neighbors",
          "Specifies that the only residues that do not have too many neighbors (as calculated by the given neighbor limit measure) are used for statistics."
        )
      ),
      m_NeighborLimitMeasure
      (
        new command::Parameter
        (
          "neighbor_limit_measure",
          "\tThe method that should be used to calculate the neighbors of residues.",
          command::ParameterCheckAllowed
          (
            storage::Vector< std::string>::Create
            (
              "NeighborCount",
              "NeighborVector"
            )
          ), "NeighborCount"
        )
      ),
      m_NeighborLimit
      (
        new command::Parameter
        (
          "neighbor_limit",
          "\tThe value which is the maximum (not inclusive) desired allowable neighbor limit.",
          "0.0"
        )
      ),
      m_NeighborLimitSeqExcl
      (
        new command::Parameter
        (
          "neighbor_limit_seq_excl",
          "\tThe integer amount of sequence exclusion that should be used to calculate neighbor measure.",
          command::ParameterCheckRanged< size_t>( 0, 9999),
          "0"
        )
      ),
      m_CalculateOnlyAccessibilities
      (
        new command::Parameter
        (
          "calculate_only_accessibilities",
          "\tOnly accessibilities will be calculated and output using BCL Message with a level of Debug..",
          command::ParameterCheckAllowed
          (
            storage::Vector< std::string>::Create
            (
              "true",
              "false"
            )
          ), "false"
        )
      ),
      m_HistogramOutputFile
      (
        new command::FlagStatic
        (
          "output_file",
          "Path and name of the output file which will hold the results.",
          command::Parameter
          (
            "string",
            "\tpath and name of the outputfile",
            "spin_label_statistics.txt"
          )
        )
      ),
      m_HistogramSpecifications
      (
        new command::FlagStatic
        (
          "histogram_specifications",
          "The specifications for the histogram."
        )
      ),
      m_LowerLimit
      (
        new command::Parameter
        (
          "lower_limit",
          "\tThe value at which the first bin of the histogram should start.",
          "0"
        )
      ),
      m_BinSize
      (
        new command::Parameter
        (
          "bin_size",
          "\tThe size each bin of the histogram should be (i.e. the width each bin should encompass).",
          "1.0"
        )
      ),
      m_NumberOfBins
      (
        new command::Parameter
        (
          "number_of_bins",
          "\tThe number of bins the histogram should contain.",
          "20"
        )
      ),
      m_SpinLabelLength
      (
        new command::FlagStatic
        (
          "sl_length",
          "\tlength of the spin label",
          command::Parameter
          (
            "length",
            "\tlength in angstrom of the spin label",
            "6.0"
          )
        )
      ),
      m_MinCACBSLAngle
      (
        new command::FlagStatic
        (
          "min_cacbsl_angle",
          "\tminimum angle CA->CB->SL can be",
          command::Parameter
          (
            "angle in radians",
            "\tmin angle for CA->CB->SL",
            util::Format()( 3.0 / 4.0 * math::g_Pi)
          )
        )
      ),
      m_MaxCACBSLAngle
      (
        new command::FlagStatic
        (
          "max_cacbsl_angle",
          "\tmaximum angle CA->CB->SL can be",
          command::Parameter
          (
            "angle in radians",
            "\tmax angle for CA->CB->SL",
            util::Format()( math::g_Pi)
          )
        )
      ),
      m_MinSLDihedralAngle
      (
        new command::FlagStatic
        (
          "min_sl_dihedral_angle",
          "\tminimum angle sl dihedral can be",
          command::Parameter
          (
            "angle in radians",
            "\tminimum angle sl dihedral can be",
            "0.0"
          )
        )
      ),
      m_MaxSLDihedralAngle
      (
        new command::FlagStatic
        (
          "max_sl_dihedral_angle",
          "\tmaximum angle sl dihedral can be",
          command::Parameter
          (
            "angle in radians",
            "\tmaximum angle sl dihedral can be",
            util::Format()( 2.0 * math::g_Pi)
          )
        )
      ),
      m_CalculateRadiusOfGyration
      (
        new command::FlagDynamic
        (
          "calc_radius_gyration",
          "\tOnly calculate the radius of gyration. No other statistics will be done. Will be put into \"output_file\".",
          command::Parameter
          (
            "backboneatom",
            "any backbone atoms from the list",
            command::ParameterCheckAllowed
            (
              biol::GetAtomTypes().GetBackBoneAtomNames()
            )
          ),
          0, 5
        )
      ),
      m_CalculateSL_CB_SL_CB
      (
        new command::FlagStatic
        (
          "Dsl-Dsb-Dsl-Dcb",
          "Do statistics involving calculations of (Dsl-Dcb) - (Dsl-Dcb)"
        )
      ),
      m_StartOffset
      (
        new command::Parameter
        (
          "start_offset",
          "\tthe number of residues to go backwards or forwards to reach the starting point",
          "0"
        )
      ),
      m_NumberOFSteps
      (
        new command::Parameter
        (
          "steps",
          "\tthe number of residues/steps to go forward",
          "0"
        )
      ),
      m_CaCbSlAngleDrift
      (
        new command::Parameter
        (
          "ca_cb_sl_angle_drift",
          "\t In radians, the amount that the walking SL CA->CB->SL angle can differ (+ and -) from the initial SL. Pi = no relation",
          util::Format()( math::g_Pi)
        )
      ),
      m_SlDihedralAngleDrift
      (
        new command::Parameter
        (
          "sl_dihedral_angle_drift",
          "\tIn radians, the amount that the walking SL dihedral angle (ProtCent->CA->CB->SL) can differ (+ and -) from the initial SL. Pi = no relation",
          util::Format()( math::g_Pi)
        )
      ),
      m_RecordDslVsDCB
      (
        new command::Parameter
        (
          "record_dsl_vs_dcb",
          "\tWhether or not Dslinitial-Dslwalking and Dcbinitial-Dcbwalking values should be stored and output into a file with name denoted by \"_dsl_vs_dcb_correlations.txt\" at the end",
          command::ParameterCheckAllowed
          (
            storage::Vector< std::string>::Create
            (
              "true",
              "false"
            )
          ),
          "false"
        )
      ),
      m_CalculateSL_CB
      (
        new command::FlagStatic
        (
          "Dsl-Dsb",
          "Do statistics involving calculations of (Dsl-Dcb)"
        )
      ),
      m_LowerLimitX
      (
        new command::Parameter
        (
          "lower_limit_x",
          "\tlowerlimit of the histogram2D in the X direction",
          "0"
        )
      ),
      m_LowerLimitY
      (
        new command::Parameter
        (
          "lower_limit_y",
          "\tlowerlimit of the histogram2D in the Y direction",
          "0"
        )
      ),
      m_BinSizeX
      (
        new command::Parameter
        (
          "bin_size_x",
          "\tsize of the bins of the histogram2D in the X direction",
          "0"
        )
      ),
      m_BinSizeY
      (
        new command::Parameter
        (
          "bin_size_y",
          "\tsize of the bins of the histogram2D in the Y direction",
          "0"
        )
      ),
      m_NumberOfBinsX
      (
        new command::Parameter
        (
          "number_bins_x",
          "\tthe number of bins in the histogram2D in the X direction",
          "0"
        )
      ),
      m_NumberOfBinsY
      (
        new command::Parameter
        (
          "number_bins_y",
          "\tthe number of bins in the histogram2D in the Y direction",
          "0"
        )
      ),
      m_CalculateExperimentalSL_CB
      (
        new command::FlagStatic
        (
          "Dsl-Dcb_experimental",
          "Do statistics involving calculations of experimentally determined (Dsl-Dcb). TO BE USED ONLY IN CONJUNCTION WITH \"-Dsl-Dsb\" FLAG"
        )
      ),
      m_RestraintFile
      (
        new command::Parameter
        (
          "restraint_filename",
          "\tthe file with the experimental data",
          "default.bcl_restraint"
        )
      ),
      m_SSETypes_Constant
      (
        new command::FlagDynamic
        (
          "sse_types_constant",
          "\tflag for specifying the types of sses for Dsl-Dcb-Dsl-DcbSSESpecific (only!) calculations for the spin label that is not stepped",
          command::Parameter
          (
            "sse_type",
            "type of sse from the list",
            command::ParameterCheckEnumerate< biol::SSTypes>()
          ), 0, biol::GetSSTypes().COIL.GetIndex() + 1
        )
      ),
      m_SSETypes_Step
      (
        new command::FlagDynamic
        (
          "sse_types_step",
          "\tflag for specifying the types of sses for Dsl-Dcb-Dsl-DcbSSESpecific (only!) calculations for the spin label that is stepped",
          command::Parameter
          (
            "sse_type",
            "type of sse from the list",
            command::ParameterCheckEnumerate< biol::SSTypes>()
          ), 0, biol::GetSSTypes().COIL.GetIndex() + 1
        )
      ),
      m_CalculateAccessibilityStatistics
      (
        new command::FlagStatic
        (
          "accessibility_statistics",
          "Do statistics involving accessibility calculation SLaccessibility-CBaccessibility"
        )
      ),
      m_AccessibilityMeasure
      (
        new command::Parameter
        (
          "accessibility_measure",
          "\tThe method that should be used to calculate the accessibility.",
          command::ParameterCheckAllowed
          (
            storage::Vector< std::string>::Create
            (
              "NeighborCount",
              "NeighborVector"
            )
          ),
          "NeighborCount"
        )
      ),
      m_AccessibilityMeasureSeqExcl
      (
        new command::Parameter
        (
          "accessibility_measure_seq_excl",
          "\tThe integer amount of sequence exclusion that should be used to calculate the accessibility statistics.",
          "0"
        )
      ),
      m_CalculateDihedralAngle
      (
        new command::FlagStatic
        (
          "calc_dihedral",
          "Only calculate a dihedral angle between four points A->B-x->C->D given four sets of coordinates"
        )
      ),
      m_CalculateProjectionAngle
      (
        new command::FlagStatic
        (
          "calc_projection_angle",
          "Only calculate a projection angle between four points B->AC->D given four sets of coordinates"
        )
      ),
      m_XCoordinateA
      (
        new command::Parameter
        (
          "x_coord_point_a",
          "\tThe x-coordinate of point A",
          "0.0"
        )
      ),
      m_YCoordinateA
      (
        new command::Parameter
        (
          "y_coord_point_a",
          "\tThe y-coordinate of point A",
          "0.0"
        )
      ),
      m_ZCoordinateA
      (
        new command::Parameter
        (
          "z_coord_point_a",
          "\tThe z-coordinate of point A",
          "0.0"
        )
      ),
      m_XCoordinateB
      (
        new command::Parameter
        (
          "x_coord_point_b",
          "\tThe x-coordinate of point B",
          "0.0"
        )
      ),
      m_YCoordinateB
      (
        new command::Parameter
        (
          "y_coord_point_b",
          "\tThe y-coordinate of point B",
          "0.0"
        )
      ),
      m_ZCoordinateB
      (
        new command::Parameter
        (
          "z_coord_point_b",
          "\tThe z-coordinate of point B",
          "0.0"
        )
      ),
      m_XCoordinateC
      (
        new command::Parameter
        (
          "x_coord_point_c",
          "\tThe x-coordinate of point C",
          "0.0"
        )
      ),
      m_YCoordinateC
      (
        new command::Parameter
        (
          "y_coord_point_c",
          "\tThe y-coordinate of point C",
          "0.0"
        )
      ),
      m_ZCoordinateC
      (
        new command::Parameter
        (
          "z_coord_point_c",
          "\tThe z-coordinate of point C",
          "0.0"
        )
      ),
      m_XCoordinateD
      (
        new command::Parameter
        (
          "x_coord_point_d",
          "\tThe x-coordinate of point D",
          "0.0"
        )
      ),
      m_YCoordinateD
      (
        new command::Parameter
        (
          "y_coord_point_d",
          "\tThe y-coordinate of point D",
          "0.0"
        )
      ),
      m_ZCoordinateD
      (
        new command::Parameter
        (
          "z_coord_point_d",
          "\tThe z-coordinate of point D",
          "0.0"
        )
      ),
      m_CreateHistogram2D
      (
        new command::FlagStatic
        (
          "create_histogram_2d",
          "Only read in value pairs and create a histogram 2d out of them which will be written to file along with gnuplot file"
        )
      ),
      m_Histogram2DValuePairFilename
      (
        new command::Parameter
        (
          "value_pair_filename",
          "\tThe file containing the value pairs which will be used to create the histogram 2d. The columns containing the x y values are specified by parameters below",
          "values.txt"
        )
      ),
      m_Histogram2DOutputFilename
      (
        new command::Parameter
        (
          "histogram_2d_output_filename",
          "\tThe filename to which the histogram 2d should be written to",
          "histogram_2d.txt"
        )
      ),
      m_Histogram2DGnuplotOutputFilename
      (
        new command::Parameter
        (
          "gnuplot_output_filename",
          "\tThe filename to which the histogram 2d gnuplot file should be written to",
          "histogram_2d_gnuplot.txt"
        )
      ),
      m_Histogram2DMinX
      (
        new command::Parameter
        (
          "min_value_x",
          "\tThe starting minimum x value",
          "0"
        )
      ),
      m_Histogram2DMinY
      (
        new command::Parameter
        (
          "min_value_y",
          "\tThe starting minimum y value",
          "0"
        )
      ),
      m_Histogram2DBinSizeX
      (
        new command::Parameter
        (
          "bin_size_x",
          "\tThe size the x bins should have",
          "1"
        )
      ),
      m_Histogram2DBinSizeY
      (
        new command::Parameter
        (
          "bin_size_y",
          "\tThe size the y bins should have",
          "1"
        )
      ),
      m_Histogram2DNumBinsX
      (
        new command::Parameter
        (
          "num_bins_x",
          "\tThe number of bins in the x direction",
          "1"
        )
      ),
      m_Histogram2DNumBinsY
      (
        new command::Parameter
        (
          "num_bins_y",
          "\tThe number of bins in the y direction",
          "1"
        )
      ),
      m_Histogram2DGnuplotTitle
      (
        new command::Parameter
        (
          "gnuplot_title",
          "\tThe title that the gnuplot should have",
          "Title"
        )
      ),
      m_Histogram2DValueColumnX
      (
        new command::Parameter
        (
          "x_value_column",
          "\tThe column in the value input file which has the x values. First column = 0, second column = 1, etc...",
          "0"
        )
      ),
      m_Histogram2DValueColumnY
      (
        new command::Parameter
        (
          "y_value_column",
          "\tThe column in the value input file which has the y values. First column = 0, second column = 1, etc...",
          "1"
        )
      ),
      m_Histogram2DNormalize
      (
        new command::Parameter
        (
          "normalize_histogram",
          "\tdetermines if the histogram is normalized or not (i.e. the sum of all counts equals 1 or not)",
          "false"
        )
      ),
      m_CreateHistogram
      (
        new command::FlagStatic
        (
          "create_histogram",
          "Only read in values and create a histogram out of them which will be written to file along with gnuplot file"
        )
      ),
      m_HistogramValueFilename
      (
        new command::Parameter
        (
          "value_filename",
          "\tThe file containing the values which will be used to create the histogram. The column containing the values are specified by parameter below",
          "values.txt"
        )
      ),
      m_HistogramOutputFilename
      (
        new command::Parameter
        (
          "histogram_output_filename",
          "\tThe filename to which the histogram should be written to",
          "histogram_.txt"
        )
      ),
      m_HistogramGnuplotOutputFilename
      (
        new command::Parameter
        (
          "gnuplot_output_filename",
          "\tThe filename to which the histogram gnuplot file should be written to",
          "histogram_gnuplot.txt"
        )
      ),
      m_HistogramMin
      (
        new command::Parameter
        (
          "min_value",
          "\tThe starting minimum value of the histogram",
          "0"
        )
      ),
      m_HistogramBinSize
      (
        new command::Parameter
        (
          "bin_size",
          "\tThe size the bins should have",
          "1"
        )
      ),
      m_HistogramNumBins
      (
        new command::Parameter
        (
          "num_bins",
          "\tThe number of bins",
          "1"
        )
      ),
      m_HistogramGnuplotTitle
      (
        new command::Parameter
        (
          "gnuplot_title",
          "\tThe title that the gnuplot should have",
          "Title"
        )
      ),
      m_HistogramValueColumn
      (
        new command::Parameter
        (
          "value_column",
          "\tThe column in the value input file which has the values. First column = 0, second column = 1, etc...",
          "0"
        )
      ),
      m_HistogramWeightColumn
      (
        new command::Parameter
        (
          "weight_column",
          "\tThe column in the input file which has the weights of the values. First column = 0, second column = 1, etc...",
          util::Format()( util::GetUndefinedSize_t())
        )
      ),
      m_HistogramNormalize
      (
        new command::Parameter
        (
          "normalize_histogram",
          "\t\"true\" or \"false\"determines if the histogram is normalized or not (i.e. the sum of all counts equals 1 or not)",
          "false"
        )
      ),
      m_NoSLClashing
      (
        new command::FlagStatic
        (
          "no_clashing",
          "Don't allow the position of the effective spin label to clash with other atoms of the protein"
        )
      ),
      m_ClashThreshold
      (
        new command::Parameter
        (
          "clash_threshold",
          "\tthe minimum radius for the SL effective that does not clash",
          "2.0"
        )
      ),
      m_ClashMaxTrials
      (
        new command::Parameter
        (
          "max_clash_trials",
          "\tthe number maximum number of trials that will be attempted to find a non clashing position for the SL effective",
          "10"
        )
      ),
      m_CalculateEPRDistanceAgreement
      (
        new command::FlagStatic
        (
          "calculate_restraint_agreement",
          "Only calculate the agreement of the list of protein models with a restraint file"
        )
      ),
      m_DistanceRestraintFile
      (
        new command::Parameter
        (
          "restraint_file",
          "\tthe restraint file that the models will be compared against",
          "restraint_file.cst"
        )
      ),
      m_DistanceRestraintScoreOutputFile
      (
        new command::Parameter
        (
          "score_output_file",
          "\tthe file which will contain the scores for the protein models",
          "output.sc"
        )
      ),
      m_EPRDistanceAgreementPDBList
      (
        new command::Parameter
        (
          "list_of_pdb_to_score",
          "\tthe file which has the list of pdbs that will be scored against the distance restraints",
          "models_to_score.ls"
        )
      ),
      m_CalculateDistanceDifferences
      (
        new command::Parameter
        (
          "calc_dist_diff",
          "\tif \"true\" will fill table with model Dcb - epr Dsl values. if \"false\" will use score",
          "false"
        )
      ),
      m_CalcSL_HStatistics
      (
        new command::FlagStatic
        (
          "calc_sl_h_cb_h_stats",
          "\tcalculate statistics for SL->H - CB->H for the protein list given"
        )
      )
    {
      // attach parameters to flags
            m_LimitNeighbors->PushBack( m_NeighborLimitMeasure);
      m_LimitNeighbors->PushBack( m_NeighborLimit);
      m_LimitNeighbors->PushBack( m_NeighborLimitSeqExcl);
      m_LimitNeighbors->PushBack( m_CalculateOnlyAccessibilities);
      m_HistogramSpecifications->PushBack( m_LowerLimit);
      m_HistogramSpecifications->PushBack( m_BinSize);
      m_HistogramSpecifications->PushBack( m_NumberOfBins);
      m_CalculateSL_CB_SL_CB->PushBack( m_StartOffset);
      m_CalculateSL_CB_SL_CB->PushBack( m_NumberOFSteps);
      m_CalculateSL_CB_SL_CB->PushBack( m_CaCbSlAngleDrift);
      m_CalculateSL_CB_SL_CB->PushBack( m_SlDihedralAngleDrift);
      m_CalculateSL_CB_SL_CB->PushBack( m_RecordDslVsDCB);
      m_CalculateSL_CB->PushBack( m_LowerLimitX);
      m_CalculateSL_CB->PushBack( m_LowerLimitY);
      m_CalculateSL_CB->PushBack( m_BinSizeX);
      m_CalculateSL_CB->PushBack( m_BinSizeY);
      m_CalculateSL_CB->PushBack( m_NumberOfBinsX);
      m_CalculateSL_CB->PushBack( m_NumberOfBinsY);
      m_CalculateExperimentalSL_CB->PushBack( m_RestraintFile);
      m_CalculateAccessibilityStatistics->PushBack( m_AccessibilityMeasure);
      m_CalculateAccessibilityStatistics->PushBack( m_AccessibilityMeasureSeqExcl);
      m_CalculateDihedralAngle->PushBack( m_XCoordinateA);
      m_CalculateDihedralAngle->PushBack( m_YCoordinateA);
      m_CalculateDihedralAngle->PushBack( m_ZCoordinateA);
      m_CalculateDihedralAngle->PushBack( m_XCoordinateB);
      m_CalculateDihedralAngle->PushBack( m_YCoordinateB);
      m_CalculateDihedralAngle->PushBack( m_ZCoordinateB);
      m_CalculateDihedralAngle->PushBack( m_XCoordinateC);
      m_CalculateDihedralAngle->PushBack( m_YCoordinateC);
      m_CalculateDihedralAngle->PushBack( m_ZCoordinateC);
      m_CalculateDihedralAngle->PushBack( m_XCoordinateD);
      m_CalculateDihedralAngle->PushBack( m_YCoordinateD);
      m_CalculateDihedralAngle->PushBack( m_ZCoordinateD);
      m_CalculateProjectionAngle->PushBack( m_XCoordinateA);
      m_CalculateProjectionAngle->PushBack( m_YCoordinateA);
      m_CalculateProjectionAngle->PushBack( m_ZCoordinateA);
      m_CalculateProjectionAngle->PushBack( m_XCoordinateB);
      m_CalculateProjectionAngle->PushBack( m_YCoordinateB);
      m_CalculateProjectionAngle->PushBack( m_ZCoordinateB);
      m_CalculateProjectionAngle->PushBack( m_XCoordinateC);
      m_CalculateProjectionAngle->PushBack( m_YCoordinateC);
      m_CalculateProjectionAngle->PushBack( m_ZCoordinateC);
      m_CalculateProjectionAngle->PushBack( m_XCoordinateD);
      m_CalculateProjectionAngle->PushBack( m_YCoordinateD);
      m_CalculateProjectionAngle->PushBack( m_ZCoordinateD);
      m_CreateHistogram2D->PushBack( m_Histogram2DValuePairFilename);
      m_CreateHistogram2D->PushBack( m_Histogram2DOutputFilename);
      m_CreateHistogram2D->PushBack( m_Histogram2DGnuplotOutputFilename);
      m_CreateHistogram2D->PushBack( m_Histogram2DMinX);
      m_CreateHistogram2D->PushBack( m_Histogram2DMinY);
      m_CreateHistogram2D->PushBack( m_Histogram2DBinSizeX);
      m_CreateHistogram2D->PushBack( m_Histogram2DBinSizeY);
      m_CreateHistogram2D->PushBack( m_Histogram2DNumBinsX);
      m_CreateHistogram2D->PushBack( m_Histogram2DNumBinsY);
      m_CreateHistogram2D->PushBack( m_Histogram2DGnuplotTitle);
      m_CreateHistogram2D->PushBack( m_Histogram2DValueColumnX);
      m_CreateHistogram2D->PushBack( m_Histogram2DValueColumnY);
      m_CreateHistogram2D->PushBack( m_Histogram2DNormalize);
      m_CreateHistogram->PushBack( m_HistogramValueFilename);
      m_CreateHistogram->PushBack( m_HistogramOutputFilename);
      m_CreateHistogram->PushBack( m_HistogramGnuplotOutputFilename);
      m_CreateHistogram->PushBack( m_HistogramMin);
      m_CreateHistogram->PushBack( m_HistogramBinSize);
      m_CreateHistogram->PushBack( m_HistogramNumBins);
      m_CreateHistogram->PushBack( m_HistogramGnuplotTitle);
      m_CreateHistogram->PushBack( m_HistogramValueColumn);
      m_CreateHistogram->PushBack( m_HistogramWeightColumn);
      m_CreateHistogram->PushBack( m_HistogramNormalize);
      m_NoSLClashing->PushBack( m_ClashThreshold);
      m_NoSLClashing->PushBack( m_ClashMaxTrials);
      m_CalculateEPRDistanceAgreement->PushBack( m_DistanceRestraintFile);
      m_CalculateEPRDistanceAgreement->PushBack( m_DistanceRestraintScoreOutputFile);
      m_CalculateEPRDistanceAgreement->PushBack( m_EPRDistanceAgreementPDBList);
      m_CalculateEPRDistanceAgreement->PushBack( m_CalculateDistanceDifferences);
    }

    const ApplicationType StatisticSpinLabel::StatisticSpinLabel_Instance
    (
      GetAppGroups().AddAppToGroup( new StatisticSpinLabel(), GetAppGroups().e_InternalBiol)
    );

  } // namespace app
} // namespace bcl

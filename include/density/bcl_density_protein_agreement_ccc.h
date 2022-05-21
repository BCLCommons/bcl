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

#ifndef BCL_DENSITY_PROTEIN_AGREEMENT_CCC_H_
#define BCL_DENSITY_PROTEIN_AGREEMENT_CCC_H_

// include the namespace header
#include "bcl_density.h"

// include other forward headers - sorted alphabetically

// includes from bcl - sorted alphabetically
#include "bcl_density_protein_agreement_interface.h"
#include "bcl_density_simulate_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinAgreementCCC
    //! @brief This class is a Deviation class
    //! @details It is derived from FunctionInterface and you can pass a densitymap to the constructor.
    //! it will the calculate the standarddeviation to a given ARGUMENT densitymap
    //!
    //! @see @link example_density_protein_agreement_ccc.cpp @endlink
    //! @author woetzen
    //! @date 08.04.2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinAgreementCCC :
      public ProteinAgreementInterface
    {
    private:

    //////////
    // data //
    //////////

      //! this this is the given density map to which every AGRUMENT density map is compared to
      util::SiPtr< const Map>               m_Map;

      //! simulator that create a simulated density map from a given set of atoms
      util::ShPtr< SimulateInterface> m_Simulate;

      //! cutoff below which voxels are not considered for CCC calculation
      double                          m_SimulatedContourLevelCutoff;

      //! add side chains to given protein model amino acids
      bool                            m_UseSideChains;

      //! mutliple the ccc with the number of amino acids in the protein model
      bool                            m_MultiplyWihtNumberAminoAcids;

    public:

    //////////
    // data //
    //////////

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param ADD_SIDECHAIN_ATOMS protein model will get side chains before the density simulation
      //! @param MULTIPLY_WITH_NUMBER_AAS multiply with number AAs so that the agreement scales with protein size
      ProteinAgreementCCC
      (
        const bool ADD_SIDECHAIN_ATOMS,
        const bool MULTIPLY_WITH_NUMBER_AAS
      );

      //! @brief virtual copy constructor
      //! @return pointer to a new ProteinAgreementCCC copied from this one
      ProteinAgreementCCC *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief access to the density simulator
      //! @return the density simulator used
      const util::ShPtr< SimulateInterface> &GetSimulator() const;

      //! @brief set the simulator used for density agreement, e.g. for simulating density from the protein model
      //! @param SP_SIMULATOR ShPtr to SimulatInterface
      void SetSimulator( const util::ShPtr< SimulateInterface> &SP_SIMULATOR);

      //! @brief access to the density used for agreement calculation
      //! @return SiPtr to the density
      const util::SiPtr< const Map> &GetDensity() const;

      //! @brief set the density used for agreement calculation
      //! @param SP_DENSITY SiPtr to the density map
      void SetDensityMap( const util::SiPtr< const Map> &SP_DENSITY);

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ///////////////
    // operators //
    ///////////////

      //! @brief operator() taking the a ProteinModel and returning the cross correlation coefficient
      //! @param PROTEIN_MODEL
      //! @return correlation between the member density map and a simulated density map for PROTEIN_MODEL
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // input and output //
    //////////////////////

    protected:

      //! @brief write to ostream
      //! @param OSTREAM is the output stream
      //! @param INDENT indentation
      //! @return returns the output stream
      std::ostream &Write( std::ostream &OSTREAM, const size_t INDENT) const;

      //! @brief read from istream
      //! @param ISTREAM is the input stream
      //! @return returns the input stream
      std::istream &Read( std::istream &ISTREAM);

    public:

    //////////////////////
    // helper functions //
    //////////////////////

      //! @brief calculate the cross correlation between experimental and simulated density map
      //! this ccc measure only considers voxels, if the intensity in the simulated map is above the given contour level
      //! this prevents that adjacent protein densities in experimental maps have a negative contribution to the measure
      //! @param EXPERIMENTAL_DENSITY_MAP map from experiment
      //! @param SIMULATED_DENSITY_MAP map simulated from protein structure
      //! @param CONTOUR_LEVEL_SIMULATED the minimal intensity in the simulated map voxel, to be considered for the CCC
      //! @return cross correlation coefficient
      static double CrossCorrelationCoefficient
      (
        const Map &EXPERIMENTAL_DENSITY_MAP,
        const Map &SIMULATED_DENSITY_MAP,
        const double CONTOUR_LEVEL_SIMULATED
      );

    }; //class Deviation

  } // namespace density
} // namespace bcl

#endif //BCL_DENSITY_PROTEIN_AGREEMENT_CCC_H_

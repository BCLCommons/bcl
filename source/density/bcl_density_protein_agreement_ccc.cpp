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
#include "density/bcl_density_protein_agreement_ccc.h"

// includes from bcl - sorted alphabetically
#include "assemble/bcl_assemble_protein_model.h"
#include "biol/bcl_biol_aa_side_chain_factory.h"
#include "biol/bcl_biol_atom.h"
#include "density/bcl_density_map.h"
#include "io/bcl_io_serialization.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace density
  {

  //////////
  // data //
  //////////

    //! single instance of that class
    const util::SiPtr< const util::ObjectInterface> ProteinAgreementCCC::s_Instance
    (
      GetObjectInstances().AddInstance( new ProteinAgreementCCC( false, false))
    );

  //////////////////////////////////
  // construction and destruction //
  //////////////////////////////////

    //! @brief default constructor
    //! @param ADD_SIDECHAIN_ATOMS protein model will get side chains before the density simulation
    //! @param MULTIPLY_WITH_NUMBER_AAS multiply with number AAs so that the agreement scales with protein size
    ProteinAgreementCCC::ProteinAgreementCCC
    (
      const bool ADD_SIDECHAIN_ATOMS,
      const bool MULTIPLY_WITH_NUMBER_AAS
    ) :
      m_SimulatedContourLevelCutoff( 0.0),
      m_UseSideChains( ADD_SIDECHAIN_ATOMS),
      m_MultiplyWihtNumberAminoAcids( MULTIPLY_WITH_NUMBER_AAS)
    {
    }

    //! @brief virtual copy constructor
    //! @return pointer to a new ProteinAgreementCCC copied from this one
    ProteinAgreementCCC *ProteinAgreementCCC::Clone() const
    {
      return new ProteinAgreementCCC( *this);
    }

  /////////////////
  // data access //
  /////////////////

    //! @brief returns class name
    //! @return the class name as const ref std::string
    const std::string &ProteinAgreementCCC::GetClassIdentifier() const
    {
      return GetStaticClassName( *this);
    }

    //! @brief gets the scheme for this function class
    //! @return the scheme for this function class
    const std::string &ProteinAgreementCCC::GetScheme() const
    {
      static const std::string s_scheme( "CCC");
      static const std::string s_scheme_sclaed( "CCCScaled");
      return m_MultiplyWihtNumberAminoAcids ? s_scheme_sclaed : s_scheme;
    }

    //! @brief access to the density simulator
    //! @return the density simulator used
    const util::ShPtr< SimulateInterface> &ProteinAgreementCCC::GetSimulator() const
    {
      return m_Simulate;
    }

    //! @brief set the simulator used for density agreement, e.g. for simulating density from the protein model
    //! @param SP_SIMULATOR ShPtr to SimulatInterface
    void ProteinAgreementCCC::SetSimulator( const util::ShPtr< SimulateInterface> &SP_SIMULATOR)
    {
      m_Simulate = SP_SIMULATOR;
    }

    //! @brief access to the density used for agreement calculation
    //! @return SiPtr to the density
    const util::SiPtr< const Map> &ProteinAgreementCCC::GetDensity() const
    {
      return m_Map;
    }

    //! @brief set the density used for agreement calculation
    //! @param SP_DENSITY SiPtr to the density map
    void ProteinAgreementCCC::SetDensityMap( const util::SiPtr< const Map> &SP_DENSITY)
    {
      m_Map = SP_DENSITY;
    }

    //! @brief return parameters for member data that are set up from the labels
    //! @return parameters for member data that are set up from the labels
    io::Serializer ProteinAgreementCCC::GetSerializer() const
    {
      io::Serializer serializer;
      serializer.SetClassDescription
      (
        "Calculates the standard deviation between a given protein model and a density map."
      );
      // serializer.AddInitializer
      // (
      //   "map",
      //   "density map to which the protein model is compared",
      //   io::Serialization::GetAgent( &m_Map)
      // );
      // serializer.AddInitializer
      // (
      //   "simulator",
      //   "algorithm to simulate a density map fora protein to facilitate comparison for CCC calculation.",
      //   io::Serialization::GetAgent( &m_Simulate)
      // );
      serializer.AddInitializer
      (
        "contour cutoff",
        "cutoff below which voxels are not considered for CCC calculation.",
        io::Serialization::GetAgent( &m_SimulatedContourLevelCutoff)
      );
      serializer.AddInitializer
      (
        "use side chains",
        "whether to add side chains to the protein model prior to CCC calculation.",
        io::Serialization::GetAgent( &m_UseSideChains)
      );
      serializer.AddInitializer
      (
        "multiply ccc",
        "multiply CCC with the number of amino acids in the protein model.",
        io::Serialization::GetAgent( &m_MultiplyWihtNumberAminoAcids)
      );

      return serializer;
    }

  ///////////////
  // operators //
  ///////////////

    //! @brief operator() taking the a ProteinModel and returning the cross correlation coefficient
    //! @param PROTEIN_MODEL
    //! @return correlation between the member density map and a simulated density map for PROTEIN_MODEL
    double ProteinAgreementCCC::operator()( const assemble::ProteinModel &PROTEIN_MODEL) const
    {
      // density map
      Map simulated_map;

      if( m_UseSideChains)
      {
        // generate new protein model with side chains from PROTEINMODEL
        const util::ShPtr< assemble::ProteinModel> protein_model_with_side_chains
        (
          biol::AASideChainFactory( false, true).ProteinModelWithSideChains( PROTEIN_MODEL)
        );

        // simulate map
        simulated_map = m_Simulate->operator ()( protein_model_with_side_chains->GetAtoms());
      }
      else
      {
        // use protein as it is
        simulated_map = m_Simulate->operator ()( PROTEIN_MODEL.GetAtoms());
      }

      // cross correlation coefficient between simulated and actual map
      double ccc( CrossCorrelationCoefficient( *m_Map, simulated_map, m_SimulatedContourLevelCutoff));

      // scale to protein size
      if( m_MultiplyWihtNumberAminoAcids)
      {
        ccc *= PROTEIN_MODEL.GetNumberAAs();
      }

      // return correlation
      return -ccc;
    }

  //////////////////////
  // input and output //
  //////////////////////

    //! @brief write to ostream
    //! @param OSTREAM is the output stream
    //! @param INDENT indentation
    //! @return returns the output stream
    std::ostream &ProteinAgreementCCC::Write( std::ostream &OSTREAM, const size_t INDENT) const
    {
      // write member
      io::Serialize::Write( m_Map                         , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_Simulate                    , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_SimulatedContourLevelCutoff , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_UseSideChains               , OSTREAM, INDENT) << '\n';
      io::Serialize::Write( m_MultiplyWihtNumberAminoAcids, OSTREAM, INDENT);

      // end
      return OSTREAM;
    }

    //! @brief read from istream
    //! @param ISTREAM is the input stream
    //! @return returns the input stream
    std::istream &ProteinAgreementCCC::Read( std::istream &ISTREAM)
    {
      // write member
      io::Serialize::Read( m_Map                         , ISTREAM);
      io::Serialize::Read( m_Simulate                    , ISTREAM);
      io::Serialize::Read( m_SimulatedContourLevelCutoff , ISTREAM);
      io::Serialize::Read( m_UseSideChains               , ISTREAM);
      io::Serialize::Read( m_MultiplyWihtNumberAminoAcids, ISTREAM);

      // end
      return ISTREAM;
    }

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
    double ProteinAgreementCCC::CrossCorrelationCoefficient
    (
      const Map &EXPERIMENTAL_DENSITY_MAP,
      const Map &SIMULATED_DENSITY_MAP,
      const double CONTOUR_LEVEL_SIMULATED
    )
    {
      // create common sub tensor
      const storage::VectorND< 2, math::Tensor< double> > exp_sim_sub_tensor( EXPERIMENTAL_DENSITY_MAP.CommonSubTensor( SIMULATED_DENSITY_MAP));

      // number of voxels that are above the CONTOUR_LEVEL in the experimental (this) density map
      size_t count_voxel( 0);
      double sum_sim( 0);
      double sum_exp( 0);
      double sum_sim2( 0);
      double sum_exp2( 0);
      double sum_exp_sim( 0);

      // iterate over experimental and simulated sub tensor
      for
      (
        const double
          *exp( exp_sim_sub_tensor.First().Begin()), *exp_end( exp_sim_sub_tensor.First().End()),
          *sim( exp_sim_sub_tensor.Second().Begin());
        exp != exp_end;
        ++exp, ++sim
      )
      {
        const double sim_int( *sim);

        // ignore sim and exp intensities below sim contour level
        if( sim_int > CONTOUR_LEVEL_SIMULATED)
        {
          const double exp_int( *exp);
          ++count_voxel;
          sum_exp += exp_int;
          sum_sim += sim_int;
          sum_exp_sim += exp_int * sim_int;
          sum_exp2 += math::Sqr( exp_int);
          sum_sim2 += math::Sqr( sim_int);
        }
      }

      // calculate actual correlation
      double correlation( count_voxel * sum_exp_sim - sum_exp * sum_sim);
      correlation /= math::Sqrt( count_voxel * sum_exp2 - math::Sqr( sum_exp)) * math::Sqrt( count_voxel * sum_sim2 - math::Sqr( sum_sim));

      // end
      return correlation;
    }

  } // namespace density
} // namespace bcl

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

#ifndef BCL_OPENCL_PROTEIN_AGREEMENT_CCC_H_
#define BCL_OPENCL_PROTEIN_AGREEMENT_CCC_H_

// include the namespace header
#include "bcl_opencl.h"

// include other forward headers - sorted alphabetically
#include "density/bcl_density.fwd.hh"
#include "linal/bcl_linal.fwd.hh"
#include "math/bcl_math.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_opencl_command_queue.h"
#include "bcl_opencl_kernel_source_interface.h"
#include "bcl_opencl_vector.h"
#include "density/bcl_density_protein_agreement_interface.h"
#include "density/bcl_density_simulate_interface.h"
#include "util/bcl_util_sh_ptr.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace opencl
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class ProteinAgreementCCC
    //! @brief This class is a Deviation class
    //! @details It is derived from FunctionInterface and you can pass a density map to the constructor.
    //! it will the calculate the standard deviation to a given ARGUMENT density map
    //!
    //! @see @link example_opencl_protein_agreement_ccc.cpp @endlink
    //! @author woetzen
    //! @date 08.04.2007
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API ProteinAgreementCCC :
      public density::ProteinAgreementInterface
    {
    private:

    //////////
    // data //
    //////////

      //! this this is the given density map to which every AGRUMENT density map is compared to
      util::SiPtr< const density::Map>         m_Map;

      //! simulator that create a simulated density map from a given set of atoms
      util::ShPtr< density::SimulateInterface> m_Simulate;

      //! cutoff below which voxels are not considered for CCC calculation
      double                          m_SimulatedContourLevelCutoff;

      //! add sidechains to given protein model
      bool m_UseSideChains;

      CommandQueue   m_CommandQueue; //!< the command queue
      cl::Program    m_Program;      //!< the opencl program

    public:

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief default constructor
      //! @param ADD_SIDECHAIN_ATOMS protein model will get side chains before the density simulation
      ProteinAgreementCCC( const bool ADD_SIDECHAIN_ATOMS = false);

      //! @brief virtual copy constructor
      //! @return pointer to a new ProteinAgreementCCC copied from this one
      ProteinAgreementCCC *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief access to the density simulator
      //! @return the density simulator used
      const util::ShPtr< density::SimulateInterface> &GetSimulator() const;

      //! @brief set the simulator used for density agreement, e.g. for simulating density from the protein model
      //! @param SP_SIMULATOR ShPtr to SimulatInterface
      void SetSimulator( const util::ShPtr< density::SimulateInterface> &SP_SIMULATOR);

      //! @brief access to the density used for agreement calculation
      //! @return SiPtr to the density
      const util::SiPtr< const density::Map> &GetDensity() const;

      //! @brief set the density used for agreement calculation
      //! @param SP_DENSITY SiPtr to the density map
      void SetDensityMap( const util::SiPtr< const density::Map> &SP_DENSITY);

    ////////////////
    // operations //
    ////////////////

      //! @brief is this class compatible with given command queue
      //! @param COMMAND_QUEUE the command queue this object would operate on
      //! @return bool true, if this command queue can be used, e.g. the devices support the necessary extensions
      bool IsCompatible( const CommandQueue &COMMAND_QUEUE) const;

      //! @brief initialize this class
      //! @brief COMMAND_QUEUE queue to use
      //! @return bool true is initialization was successful, false otherwise (when program could not be compiled ...)
      bool Initialize( const CommandQueue &COMMAND_QUEUE);

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

      static const KernelSource &GetCCCKernel();

      //! @brief copy the density map to the device
      //! @param DENSITY_MAP the map to be copied
      //! @param NEW_DIMENSIONS new dimensions for map - used for padding
      //! @return Buffer the buffer containing the map
      Vector< double> MapToDevice
      (
        const density::Map &DENSITY_MAP,
        const storage::VectorND< 3, size_t> &NEW_DIMENSIONS
      ) const;

      //! @brief copy the tensor to the device
      //! @param TENSOR the tensor to be copied
      //! @return Buffer the buffer containing the tensor
      Vector< double> TensorToDevice( const math::Tensor< double> &TENSOR) const;

      //! @brief calculate the cross correlation between experimental and simulated density map
      //! this ccc measure only considers voxels, if the intensity in the simulated map is above the given contour level
      //! this prevents that adjacent protein densities in experimental maps have a negative contribution to the measure
      //! @param EXPERIMENTAL_DENSITY_MAP map from experiment
      //! @param SIMULATED_DENSITY_MAP map simulated from protein structure
      //! @param CONTOUR_LEVEL_SIMULATED the minimal intensity in the simulated map voxel, to be considered for the CCC
      //! @return cross correlation coefficient
      double CrossCorrelationCoefficient
      (
        const density::Map &EXPERIMENTAL_DENSITY_MAP,
        const density::Map &SIMULATED_DENSITY_MAP,
        const double CONTOUR_LEVEL_SIMULATED
      ) const;

      //! @brief calculate the cross correlation between experimental and simulated density map
      //! this ccc measure only considers voxels, if the intensity in the simulated map is above the given contour level
      //! this prevents that adjacent protein densities in experimental maps have a negative contribution to the measure
      //! @param EXPERIMENTAL_SUBDENSITY_MAP map from experiment
      //! @param SIMULATED_SUBDENSITY_MAP map simulated from protein structure
      //! @param CONTOUR_LEVEL_SIMULATED the minimal intensity in the simulated map voxel, to be considered for the CCC
      //! @return cross correlation coefficient
      double CrossCorrelationCoefficient
      (
        const math::Tensor< double> &EXPERIMENTAL_SUBDENSITY_MAP,
        const math::Tensor< double> &SIMULATED_SUBDENSITY_MAP,
        const double CONTOUR_LEVEL_SIMULATED
      ) const;

      //! @brief calculate the cross correlation between experimental and simulated density map
      //! @see CrossCorrelationCoefficient
      //! @param EXPERIMENTAL_BUFFER map from experiment
      //! @param SIMULATED_BUFFER map simulated from protein structure
      //! @param GRID_SIZE number of elements in buffer
      //! @param CONTOUR_LEVEL_SIMULATED the minimal intensity in the simulated map voxel, to be considered for the CCC
      //! @return cross correlation coefficient
      double CrossCorrelationCoefficient
      (
        const Vector< double> &EXPERIMENTAL_BUFFER,
        const Vector< double> &SIMULATED_BUFFER,
        const size_t GRID_SIZE,
        const double CONTOUR_LEVEL_SIMULATED
      ) const;

      double CrossCorrelationCoefficient
      (
        const Vector< double> &EXPERIMENTAL_BUFFER,
        const Vector< double> &SIMULATED_BUFFER,
        const linal::Vector< int> &EXP_START,
        const storage::VectorND< 3, size_t> &EXP_DIMENSIONS,
        const linal::Vector< int> &SIM_START,
        const storage::VectorND< 3, size_t> &SIM_DIMENSIONS,
        const linal::Vector< int> &EXTENSION,
        const double CONTOUR_LEVEL_SIMULATED
      ) const;

      //! @brief converts tensor< double> to tensor< float>
      math::Tensor< float> Convert( const math::Tensor< double> &TENSOR) const;

    }; // class ProteinAgreementCCC

  } // namespace opencl
} // namespace bcl

#endif //BCL_OPENCL_PROTEIN_AGREEMENT_CCC_H_

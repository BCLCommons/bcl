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

#ifndef BCL_SCORE_RADIUS_OF_GYRATION_H_
#define BCL_SCORE_RADIUS_OF_GYRATION_H_

// include the namespace header
#include "bcl_score.h"

// include other forward headers - sorted alphabetically
#include "assemble/bcl_assemble.fwd.hh"
#include "math/bcl_math.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "bcl_score_protein_model.h"
#include "biol/bcl_biol_atom_types.h"
#include "biol/bcl_biol_membrane.h"
#include "math/bcl_math_cubic_spline_damped.h"
#include "util/bcl_util_implementation.h"

// external includes - sorted alphabetically

namespace bcl
{
  namespace score
  {

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class RadiusOfGyration
    //! @brief This is a template class for scoring the radius of gyration of a set of aminoacids.\n
    //! @details It passes the CB coordinates to the function that works simply on coordinates and uses a spline to evaluate.\n
    //!
    //! @see @link example_score_radius_of_gyration.cpp @endlink
    //! @author woetzen, karakam
    //! @date 12.04.2006
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API RadiusOfGyration :
      public ProteinModel
    {

    private:

    //////////
    // data //
    //////////

      //! flag to enable normalization by sequence length
      bool m_Normalize;

      //! flag to ignore histograms -- just return raw radius of gyration
      bool m_Raw;

      //! scheme to be used in outputting
      std::string m_Scheme;

      //! path to file where the statistics and in consequence the energy potentials are read from for soluble proteins
      std::string m_HistogramFileNameSoluble;

      //! path to file where the statistics and in consequence the energy potentials are read from for membrane proteins
      std::string m_HistogramFileNameMembrane;

      //! ShPtr to the cubicspline that is used as a energy function for soluble proteins
      util::ShPtr< math::CubicSplineDamped> m_EnergyFunctionSoluble;

      //! ShPtr to the cubicspline that is used as a energy function for membrane proteins
      util::ShPtr< math::CubicSplineDamped> m_EnergyFunctionMembrane;

      //! single instance of that class
      static const util::SiPtr< const util::ObjectInterface> s_Instance;

    public:

    //////////
    // data //
    //////////

      //! @brief returns default file where the statistics and in consequence the energy potentials are read from
      //! @return default file where the statistics and in consequence the energy potentials are read from
      static const std::string &GetDefaultHistogramFilename();

      //! @brief returns default file where the statistics and in consequence the energy potentials are read from
      //! @return default file where the statistics and in consequence the energy potentials are read from
      static const std::string &GetDefaultHistogramFilenameMembrane();

    //////////////////////////////////
    // construction and destruction //
    //////////////////////////////////

      //! @brief constructor from a specified histogram file
      //! @param NORMALIZE flag to enable normalization
      //! @param RAW ignore energy function, just return raw radius of gyration
      //! @param SCHEME scheme to be used
      //! @param SOLUBLE_FILENAME filename for soluble histogram
      //! @param MEMBRANE_FILENAME filename for membrane histogram
      RadiusOfGyration
      (
        const bool NORMALIZE = false,
        const bool RAW = false,
        const std::string &SCHEME = "",
        const std::string &SOLUBLE_FILENAME = GetDefaultHistogramFilename(),
        const std::string &MEMBRANE_FILENAME = GetDefaultHistogramFilenameMembrane()
      );

      //! @brief virtual copy constructor
      RadiusOfGyration *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name
      //! @return the class name as const ref std::string
      const std::string &GetClassIdentifier() const;

      //! @brief gets the scheme for this function class
      //! @return the scheme for this function class
      const std::string &GetScheme() const;

      //! @brief returns filename of the histogram being used for soluble proteins
      //! @return filename of the histogram being used for soluble proteins
      const std::string &GetHistogramFilenameSoluble() const;

      //! @brief returns filename of the histogram being used for membrane proteins
      //! @return filename of the histogram being used for membrane proteins
      const std::string &GetHistogramFilenameMembrane() const;

      //! @brief returns the soluble energy function
      //! @return soluble energy function
      const math::CubicSplineDamped &GetEnergyFunctionSoluble() const;

      //! @brief returns the membrane energy function
      //! @return membrane energy function
      const math::CubicSplineDamped &GetEnergyFunctionMembrane() const;

      //! @brief get a more readable score scheme
      //! @return a more readable score scheme
      const std::string &GetReadableScheme() const;

      //! @brief get score type
      //! @return score type
      ProteinModel::Type GetType() const;

      //! @brief return parameters for member data that are set up from the labels
      //! @return parameters for member data that are set up from the labels
      io::Serializer GetSerializer() const;

    ////////////////
    // operations //
    ////////////////

      //! @brief calculates the radius of gyration for CB atoms for ProteinModel
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @return the radius of gyration for CB atoms for ProteinModel
      static double SquareRadiusOfGyration( const assemble::ProteinModel &PROTEIN_MODEL);

      //! @brief calculates the radius of gyration for specified atom type for ProteinModel with a collapsing normal
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @param ATOM_TYPES Atom types that will be considered when doing the calculations
      //! @return the radius of gyration for specified atom type for ProteinModel with a collapsing normal
      static double SquareRadiusOfGyrationCollapsed
      (
        const assemble::ProteinModel &PROTEIN_MODEL,
        const storage::Set< biol::AtomType> &ATOM_TYPES
      );

      //! @brief compute the square radius of gyration for a set of coordinates and a collapsing normal
      //! @param COORDINATES coordinates to be collapsed
      //! @param MEMBRANE membrane object to collapse
      //! @return square radius of gyration collapsed
      static double SquareRadiusOfGyrationCollapsed
      (
        const util::SiPtrVector< const linal::Vector3D> &COORDINATES,
        const util::ShPtr< biol::Membrane> &MEMBRANE
      );

      //! @brief Set the members of this object from the given LABEL
      //! @param LABEL the label containing members that should be read of this class
      //! @param ERROR_STREAM stream with which to write errors
      //! @return return code indicating success or failure
      bool ReadInitializerSuccessHook( const util::ObjectDataLabel &LABEL, std::ostream &ERROR_STREAM);

    ///////////////
    // operators //
    ///////////////

      //! @brief calculate the score of radius of gyration for the given ProteinModel
      //! @param PROTEIN_MODEL ProteinModel of interest
      //! @return the score of radius of gyration for the given ProteinModel
      double operator()( const assemble::ProteinModel &PROTEIN_MODEL) const;

    //////////////////////
    // helper functions //
    //////////////////////

    public:

      //! @brief create histograms from a table generated from statistics application
      //! @param TABLE containing rows: rgyr_sqr, nr_coordinates, subunits
      //! @return map containg histograms for proteins with a single subunit and one over all with the filenames as key
      static storage::Map< std::string, math::Histogram> HistogramsFromTable( const storage::Table< double> &TABLE);

    private:

      //! @brief read individual energy functions for scoring radius of gyration
      //! @param FILENAME filename to be read in
      //! @return read in energy function
      static util::ShPtr< math::CubicSplineDamped> ReadEnergyFunction( const std::string &FILENAME);

    }; // class RadiusOfGyration

  } // namespace score
} // namespace bcl

#endif //BCL_SCORE_RADIUS_OF_GYRATION_H_

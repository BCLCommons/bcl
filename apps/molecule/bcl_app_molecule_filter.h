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

#ifndef BCL_APP_MOLECULE_FILTER_H_
#define BCL_APP_MOLECULE_FILTER_H_

// include the interface for all apps
#include "app/bcl_app_apps.h"

// include other forward headers - sorted alphabetically
#include "chemistry/bcl_chemistry.fwd.hh"
#include "graph/bcl_graph.fwd.hh"
#include "storage/bcl_storage.fwd.hh"

// includes from bcl - sorted alphabetically
#include "chemistry/bcl_chemistry_fragment_ensemble.h"
#include "command/bcl_command_command.h"
#include "descriptor/bcl_descriptor_base.h"
#include "io/bcl_io_ofstream.h"
#include "math/bcl_math_comparisons.h"
#include "storage/bcl_storage_map.h"
#include "util/bcl_util_stopwatch.h"
#include "util/bcl_util_string_replacement.h"
namespace bcl
{
  namespace app
  {
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class FilterStatistics
    //! @brief class that tracks statistics about how long each filter took and how many molecules made it through
    //! @remarks example unnecessary
    //! @author mendenjl
    //! @date May 22, 2012
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    class BCL_API FilterStatistics
    {
    public:
    ///////////
    // enums //
    ///////////

      enum Type
      {
        e_DefinedAtomTypes,              //!< defined atom types
        e_3d,                            //!< valid 3d coordinates
        e_PlanarAmideBonds,              //!< Amide bonds are planar
        e_PlanarUnbridgedAromaticRings,  //!< Aromatic rings that are expected to be planar are indeed planar
        e_Simple,                        //!< simple connectivity -> not a molecular complex
        e_ContainingFragments,           //!< at least one of the desired substructures was found
        e_Matches,                       //!< molecule is isomorphic to a molecule in the given file
        e_PropertiesConstraint,          //!< Molecule has the property given in the description
        e_ComparePropertyValues,         //!< Molecule satisfied the desired property comparison
        s_NumberFilterStatistics
      };

      //! @brief get a description for the given type
      //! @param TYPE the desired type
      //! @return a string description for that type
      static const std::string &GetTypeDescription( const Type &TYPE);

    private:

    //////////
    // data //
    //////////

      util::Stopwatch m_Timer;       //!< How much time was required by the filter
      Type            m_Type;        //!< Type of statistic
      std::string     m_Description; //!< Description of the statistic, independent of the type
      size_t          m_Count;       //!< Number of molecules that the filter caught

    public:

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! @brief default constructor
      FilterStatistics();

      //! @brief constructor from description
      FilterStatistics( const Type &TYPE, const std::string &DESCRIPTION = std::string());

      //! @brief register the start of a task
      void StartFilter();

      //! @brief stop the task, increment the count
      //! @param MATCHED did the molecule match the filter
      void StopFilter( const bool &MATCHED);

      //! @brief get the full description
      //! @return the full description
      std::string GetDescription() const
      {
        return m_Description;
      }

      //! @brief get the matched count
      //! @return the matched count
      const size_t &GetMatchedCount() const
      {
        return m_Count;
      }

      //! @brief get the timer
      //! @return the timer
      const util::Stopwatch &GetTimer() const
      {
        return m_Timer;
      }

    };

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //!
    //! @class MoleculeFilter
    //! @brief Application for filtering small molecule ensembles by some criteria
    //! @details Criteria include
    //!            * do (and/or do not) contain particular fragments (given in a file)
    //!            * do not correspond to any molecule in another ensemble (at the constitutional level)
    //!            * they contain a particular misc. property
    //!            * a particular misc property value is exactly (non)equal a string (e.g. Inhibitor)
    //!            * a particular misc property value contains a particular string (e.g. nan)
    //!            * a particular small molecule misc property value is defined and satisfies a comparison (e.g. EC50 < 0.1)
    //!            * are complexes
    //!
    //! @author mendenjl
    //! @see @link example_app_molecule_filter.cpp @endlink
    //! @date September 16, 2010
    //!
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    class BCL_API MoleculeFilter :
      public InterfaceRelease
    {
    private:

    //////////
    // data //
    //////////

      //! Where to place the molecules that satisfy all the criteria given below
      util::ShPtr< command::FlagInterface> m_OutputFilenameForMoleculesSatisfyingCriteriaFlag;

      //! Where to place the molecules that fail any of the criteria given below
      util::ShPtr< command::FlagInterface> m_OutputFilenameForMoleculesFailingCriteriaFlag;

      // Basic criteria

      //! Check that all tetrahedrally hybridized atoms with 4 neighbors has at least one neighbor w/a z-coordinate
      util::ShPtr< command::FlagInterface> m_Match3DCoordinatesFlag;

      //! Check that all atom types are defined
      util::ShPtr< command::FlagInterface> m_MatchDefinedAtomTypesFlag;

      //! Check that amide bonds are planar
      util::ShPtr< command::FlagInterface> m_MatchPlanarAmideBondsFlag;

      //! Check that aromatic rings that should be planar (e.g. not cyclophanes) are planar
      util::ShPtr< command::FlagInterface> m_MatchPlanarAromaticRingsFlag;

      //! Indicate the method by which conformers will be compared
      util::ShPtr< command::FlagInterface> m_ConformerComparerFlag;
      // Fragment critera

      //! whether to match only simple molecules, e.g. not molecular complexes
      util::ShPtr< command::FlagInterface> m_MatchSimpleMoleculesFlag;

      //! input file containing sdfs of fragments
      //! Molecules with at least 1 of the fragments satisfy the criteria
      util::ShPtr< command::FlagInterface> m_MoleculesContainingFragmentsFilenameFlag;

      //! input file containing sdfs of conformers
      //! Molecules within tolerance satisfy the criteria
      util::ShPtr< command::FlagInterface> m_MoleculesContainingConformersFilenameFlag;

      // Molecule criteria

      //! input file containing sdfs of molecules
      //! Molecules in the ensemble that are also in filename satisfy this criteria
      util::ShPtr< command::FlagInterface> m_MoleculesToFindFilenameFlag;

      // Property-based criteria

      //! Molecules with these misc propertoes satisfy this criteria
      util::ShPtr< command::FlagInterface> m_HasPropertiesFlag;

      //! 2-parameter: property-name(name of a misc property), string(the value of the property)
      //! Molecules with a specified misc property (property-name) with a value string exactly equal to string satisfy this criteria
      //! Example: -has_property_with_string MiscProperty("mGlur5 Category",1) "Inhibitor"
      util::ShPtr< command::FlagInterface> m_HasPropertyWithStringFlag;

      //! 2-parameter: property-name(name of a misc property), string(the value of the property)
      //! Molecules with a specified misc property (property-name) with a value string containing string satisfy this property
      //! Example: -has_property_with_value MiscProperty("Atom_SigmaCharge",1) nan
      //! finds all the molecules whose sigma charge (which must already be stored in the sdf file) had a nan in it
      util::ShPtr< command::FlagInterface> m_HasPropertyContainingStringFlag;

      // numeric value of property with only one value

      //! 3 parameter: property-name (name of a misc property already in the sdf file),
      //!              comparison operator to use (one of >,<,>=,<=,==,!=)
      //!              value (any numeric value)
      //! if the molecule has a misc property called property-name and its value satisfies the comparison operator when
      //! compared to value, it satisfies this property
      util::ShPtr< command::FlagInterface> m_ComparePropertyValuesFlag;

      //! Set this flag so that if *any* of the comparisons is true, the molecule will be matched (rather than the
      //! default of all)
      util::ShPtr< command::FlagInterface> m_AnyFlag;

      //! statistics for each filter
      mutable std::vector< FilterStatistics> m_FilterStatistics[ FilterStatistics::s_NumberFilterStatistics];

      //! amide bond tolerance
      mutable double m_AmideBondTolerance;

      //! property constraints
      mutable storage::List< std::pair< std::string, util::StringReplacement> > m_PropertyConstraints;

      //! property comparisons
      mutable storage::List
      <
        storage::Triplet
        <
          descriptor::CheminfoProperty,
          math::Comparisons< float>::Comparison,
          descriptor::CheminfoProperty
        >
      > m_PropertyComparisons;

      //! Substructures to look for
      mutable storage::Vector< graph::ConstGraph< size_t, size_t> > m_SearchSubstructures;

      mutable chemistry::FragmentEnsemble m_ConformersToMatch;

      //! Molecules to look for, mapped to by hash strings
      mutable storage::Map< std::string, storage::List< storage::Pair< size_t, graph::ConstGraph< size_t, size_t> > > >
        m_ExactMatches;

      //! output files
      mutable io::OFStream m_OutputMatched;
      mutable io::OFStream m_OutputUnmatched;

    ///////////////////////////////////
    // construction and destruction //
    ///////////////////////////////////

      //! default constructor
      MoleculeFilter();

      //! copy constructor; ignores everything but flags
      MoleculeFilter( const MoleculeFilter &PARENT);

    public:

      // instantiate enumerator for PrepareSmallMoleculeEnsemble class
      static const ApplicationType MoleculeFilter_Instance;

      //! @brief Clone function
      //! @return pointer to new FoldProtein
      MoleculeFilter *Clone() const;

    /////////////////
    // data access //
    /////////////////

      //! @brief returns class name of the object behind a pointer or the current object
      //! @return the class name
      const std::string &GetClassIdentifier() const;

      //! @brief get the full description
      //! @return the full description
      std::string GetDescription() const;

      //! @brief initializes the command object for that executable
      util::ShPtr< command::Command> InitializeCommand() const;

      //! @brief returns readme information
      //! @return string containing information about application
      const std::string &GetReadMe() const;

      //! @brief returns web text information
      //! @return text (html allowed but not required) that will be displayed on the website
      //! @note use confluence syntax, e.g. !my_image.png! to refer to an image under documentation/image
      const std::string &GetWebText() const;

      //! @brief the Main function
      //! @return error code - 0 for success
      int Main() const;

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

      //! @brief prepare this object to handle the flags that were passed
      void Initialize() const;

      //! @brief make a hash string from a map of size-t's to size-t's
      //! @param MAP the map to make a hashable string from
      //! @return a string containing key,value,key,value pairs
      static std::string MakeHashStringFromMap( const storage::Map< size_t, size_t> &MAP);

      //! @brief make a hash string from a graph
      //! @param GRAPH a graph
      //! This string can be used as a key to a map that holds graphs with identical hash strings, thus narrowing the
      //! number of graphs that must be searched to determine whether a new scaffold is unique
      //! @return a string that specifies something about the graph that is vertex-invariant
      //!         e.g. does not depend on the ordering of the vertices
      static std::string MakeHashStringFromGraph( const graph::ConstGraph< size_t, size_t> &GRAPH);

      //! @brief apply all the filters
      //! @return true if the molecule matched all the filters
      bool MatchesFilters( chemistry::FragmentComplete &MOLECULE) const;

    }; // MoleculeFilter

  } // namespace app
} // namespace bcl
#endif // BCL_APP_MOLECULE_FILTER_H_

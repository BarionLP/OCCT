// Created by: Eugeny MALTCHIKOV
// Copyright (c) 2018 OPEN CASCADE SAS
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

#include <BOPAlgo_RemoveFeatures.hxx>

#include <BOPAlgo_Alerts.hxx>
#include <BOPAlgo_BOP.hxx>
#include <BOPAlgo_Builder.hxx>
#include <BOPAlgo_BuilderSolid.hxx>
#include <BOPAlgo_MakerVolume.hxx>
#include <BOPAlgo_Tools.hxx>
#include <Standard_ErrorHandler.hxx>

#include <BOPTools_AlgoTools.hxx>
#include <BOPTools_AlgoTools2D.hxx>
#include <BOPTools_Set.hxx>

#include <Bnd_Box.hxx>
#include <GProp_GProps.hxx>

#include <BRep_Builder.hxx>
#include <BRep_Tool.hxx>

#include <Adaptor3d_CurveOnSurface.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <Geom2dAdaptor_Curve.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <GeomPlate_BuildPlateSurface.hxx>
#include <GeomPlate_CurveConstraint.hxx>
#include <GeomPlate_MakeApprox.hxx>
#include <GeomPlate_Surface.hxx>

#include <BRepBndLib.hxx>
#include <BRepGProp.hxx>
#include <BRepCheck_Shell.hxx>

#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>

#include <BRepLib.hxx>
#include <BRepLib_MakeFace.hxx>
#include <BRepTools.hxx>

#include <ElCLib.hxx>

#include <Geom_BezierSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Geom_Curve.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <Geom_Surface.hxx>
#include <Geom_SweptSurface.hxx>
#include <GeomLib.hxx>
#include <GeomLProp_SLProps.hxx>
#include <Geom2d_Curve.hxx>
#include <Geom2d_Line.hxx>
#include <GeomAPI_ProjectPointOnSurf.hxx>
#include <GeomConvert_SurfToAnaSurf.hxx>

#include <Geom2d_Curve.hxx>

#include <NCollection_Array1.hxx>
#include <NCollection_DataMap.hxx>
#include <NCollection_DynamicArray.hxx>

#include <ShapeUpgrade_UnifySameDomain.hxx>

#include <TopAbs_ShapeEnum.hxx>

#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>

#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Iterator.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>
#include <TopoDS_Vertex.hxx>
#include <TopLoc_Location.hxx>

#include <TopoDS_Shape.hxx>
#include <TopTools_ShapeMapHasher.hxx>
#include <NCollection_IndexedDataMap.hxx>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

//=======================================================================
// static methods declaration
//=======================================================================

static void MakeRemoved(
  const NCollection_List<TopoDS_Shape>&                                theShapes,
  BRepTools_History&                                                   theHistory,
  const NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>& theKeepShapes =
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>());

static void FindInternals(const TopoDS_Shape& theS, NCollection_List<TopoDS_Shape>& theLInt);

static void RemoveInternalWires(const NCollection_List<TopoDS_Shape>& theShapes,
                                NCollection_List<TopoDS_Shape>*       theRemoved = nullptr);

static void GetOriginalFaces(
  const TopoDS_Shape&                                                  theShape,
  const NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>& theSolids,
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>&        theFeatureFacesMap,
  const NCollection_IndexedDataMap<TopoDS_Shape,
                                   NCollection_List<TopoDS_Shape>,
                                   TopTools_ShapeMapHasher>&           theAdjFaces,
  const occ::handle<BRepTools_History>&                                theHistory,
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>&       theFacesToBeKept,
  NCollection_List<TopoDS_Shape>&                                      theInternalShapes,
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>&              theFacesToCheckOri,
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>&       theSolidsToRebuild,
  NCollection_List<TopoDS_Shape>&                                      theSharedFaces,
  NCollection_List<TopoDS_Shape>&                                      theUnTouchedSolids);

static TopoDS_Shape RebuildSolidLocally(
  const TopoDS_Shape&                                           theSolid,
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theFeatureFacesMap,
  const NCollection_IndexedDataMap<TopoDS_Shape,
                                   NCollection_List<TopoDS_Shape>,
                                   TopTools_ShapeMapHasher>&    theAdjFaces,
  const occ::handle<BRepTools_History>&                         theAdjFacesHistory,
  const double                                                  theFuzzyValue,
  const bool                                                    theRunParallel,
  BOPAlgo_Builder&                                              theGF,
  const Message_ProgressRange&                                  theRange);

static void FindShape(const TopoDS_Shape& theSWhat,
                      const TopoDS_Shape& theSWhere,
                      TopoDS_Shape&       theSFound);

static void GetValidSolids(
  BOPAlgo_MakerVolume&                                          theMV,
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theFacesToCheckOri,
  const NCollection_IndexedDataMap<TopoDS_Shape,
                                   NCollection_List<TopoDS_Shape>,
                                   TopTools_ShapeMapHasher>&    theAdjFaces,
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theAnchoredFaces,
  const NCollection_List<TopoDS_Shape>&                         aSharedFaces,
  const TopoDS_Shape&                                           theOrigFaces,
  const int                                                     theNbSol,
  NCollection_List<TopoDS_Shape>&                               theLSRes,
  NCollection_List<TopoDS_Shape>&                               theRemovedShapes);

static void FindExtraShapes(
  const NCollection_IndexedDataMap<TopoDS_Shape,
                                   NCollection_List<TopoDS_Shape>,
                                   TopTools_ShapeMapHasher>&    theConnectionMap,
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theShapesToCheckOri,
  BOPAlgo_Builder&                                              theBuilder,
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>&       theShapesToAvoid,
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>*       theValidShapes = nullptr,
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>* theInternalOK  = nullptr);

static void AvoidExtraSharedFaces(NCollection_List<TopoDS_Shape>&       theLSolids,
                                  const NCollection_List<TopoDS_Shape>& theLFSharedToAvoid,
                                  BOPAlgo_Builder&                      theBuilder,
                                  NCollection_List<TopoDS_Shape>&       theExtraFaces);

static void FillSolidsHistory(
  const NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>& theSolIn,
  NCollection_List<TopoDS_Shape>&                                      theSolidsRes,
  const NCollection_IndexedDataMap<TopoDS_Shape,
                                   NCollection_List<TopoDS_Shape>,
                                   TopTools_ShapeMapHasher>&           theAdjFaces,
  BOPAlgo_Builder&                                                     theBuilder,
  BRepTools_History&                                                   theSolidsHistory);

static void TakeModified(const TopoDS_Shape&             theS,
                         BOPAlgo_Builder&                theBuilder,
                         NCollection_List<TopoDS_Shape>& theList);

static void TakeModified(const TopoDS_Shape&                                     theS,
                         BOPAlgo_Builder&                                        theBuilder,
                         NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theMap);

static void FindSolid(const TopoDS_Shape&                                        theSolIn,
                      const NCollection_List<TopoDS_Shape>&                      theSolidsRes,
                      const NCollection_IndexedDataMap<TopoDS_Shape,
                                                       NCollection_List<TopoDS_Shape>,
                                                       TopTools_ShapeMapHasher>& theAdjFaces,
                      BOPAlgo_Builder&                                           theBuilder,
                      TopoDS_Shape&                                              theSolOut);

namespace
{
//=======================================================================
// function : BOPAlgo_PIOperation
// purpose  : List of operations to be supported by the Progress Indicator
//=======================================================================
enum BOPAlgo_PIOperation
{
  PIOperation_PrepareFeatures = 0,
  PIOperation_RemoveFeatures,
  PIOperation_UpdateHistory,
  PIOperation_SimplifyResult,

  PIOperation_Last
};
} // namespace

//=================================================================================================

void BOPAlgo_RemoveFeatures::fillPIConstants(const double theWhole, BOPAlgo_PISteps& theSteps) const
{
  theSteps.SetStep(PIOperation_PrepareFeatures, 0.05 * theWhole);
  theSteps.SetStep(PIOperation_RemoveFeatures, 0.8 * theWhole);
  theSteps.SetStep(PIOperation_UpdateHistory, 0.05 * theWhole);
  theSteps.SetStep(PIOperation_SimplifyResult, 0.1 * theWhole);
}

//=======================================================================
// function: Perform
// purpose: Performs the removal of the requested faces from the input shape
//=======================================================================
void BOPAlgo_RemoveFeatures::Perform(const Message_ProgressRange& theRange)
{
  try
  {
    OCC_CATCH_SIGNALS

    if (HasHistory())
    {
      myHistory = new BRepTools_History();
    }

    // Check the input data
    CheckData();
    if (HasErrors())
    {
      return;
    }
    Message_ProgressScope aPS(theRange, "Removing features", 100);
    BOPAlgo_PISteps       aSteps(PIOperation_Last);
    analyzeProgress(100., aSteps);

    // Prepare the faces to remove.
    PrepareFeatures(aPS.Next(aSteps.GetStep(PIOperation_PrepareFeatures)));
    if (HasErrors())
    {
      return;
    }

    // Remove the features and fill the created gaps
    RemoveFeatures(aPS.Next(aSteps.GetStep(PIOperation_RemoveFeatures)));
    if (HasErrors())
    {
      return;
    }

    // Update history with the removed features
    UpdateHistory(aPS.Next(aSteps.GetStep(PIOperation_UpdateHistory)));
    if (HasErrors())
    {
      return;
    }
    // Simplify the result
    SimplifyResult(aPS.Next(aSteps.GetStep(PIOperation_SimplifyResult)));
    if (HasErrors())
    {
      return;
    }
    // Check the patches filling the gaps
    NCollection_List<TopoDS_Shape>::Iterator itP(myPatches);
    NCollection_List<double>::Iterator       itS(myPatchSizes);
    for (; itP.More() && itS.More(); itP.Next(), itS.Next())
    {
      CheckPatchTolerance(itP.Value(), itS.Value());
    }
    // Post treatment
    PostTreat();
  }
  catch (Standard_Failure const&)
  {
    AddError(new BOPAlgo_AlertRemoveFeaturesFailed());
  }
}

//=======================================================================
// function: CheckData
// purpose: Checks the input data on validity for the algorithm
//=======================================================================
void BOPAlgo_RemoveFeatures::CheckData()
{
  // Prepare the shape to work with
  myShape = myInputShape;

  // Check the type of input shape
  const TopAbs_ShapeEnum aType = myInputShape.ShapeType();

  if (aType == TopAbs_SOLID || aType == TopAbs_COMPSOLID)
  {
    return; // OK
  }

  if (aType == TopAbs_COMPOUND)
  {
    NCollection_List<TopoDS_Shape>                         aShapes;
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aMFence;
    // Extract all shapes from the compound
    BOPTools_AlgoTools::TreatCompound(myInputShape, aShapes, &aMFence);
    if (aShapes.IsEmpty())
    {
      // Add error of empty input shape
      AddError(new BOPAlgo_AlertEmptyShape(myInputShape));
      return;
    }

    // Find all solids in the list of shapes
    NCollection_List<TopoDS_Shape>           aSolids;
    NCollection_List<TopoDS_Shape>           anOtherShapes;
    NCollection_List<TopoDS_Shape>::Iterator aIt(aShapes);
    for (; aIt.More(); aIt.Next())
    {
      const TopoDS_Shape& aS = aIt.Value();
      if (aS.ShapeType() == TopAbs_SOLID || aS.ShapeType() == TopAbs_COMPSOLID)
      {
        aSolids.Append(aS);
      }
      else
      {
        anOtherShapes.Append(aS);
      }
    }

    if (aSolids.IsEmpty())
    {
      // No solids have been found for processing.
      // Add error of unsupported type of input shape
      AddError(new BOPAlgo_AlertTooFewArguments());
    }
    else if (anOtherShapes.Extent() > 0)
    {
      // Add warning of unsupported type of input shape for all
      // non-solid shapes, contained in the input shape
      for (aIt.Initialize(anOtherShapes); aIt.More(); aIt.Next())
      {
        AddWarning(new BOPAlgo_AlertUnsupportedType(aIt.Value()));
      }

      // Collect all solids into compound and overwrite the shape to rebuild
      TopoDS_Compound aCS;
      BRep_Builder().MakeCompound(aCS);
      for (aIt.Initialize(aSolids); aIt.More(); aIt.Next())
      {
        BRep_Builder().Add(aCS, aIt.Value());
      }

      myShape = aCS;

      if (HasHistory())
      {
        // Make non solid shapes removed in the history
        MakeRemoved(anOtherShapes, *myHistory.get());
      }
    }
  }
  else
  {
    // Add error of unsupported type of input shape
    AddError(new BOPAlgo_AlertTooFewArguments());
  }
}

//=======================================================================
// function: PrepareFeatures
// purpose: Prepares the features to remove
//=======================================================================
void BOPAlgo_RemoveFeatures::PrepareFeatures(const Message_ProgressRange& theRange)
{
  // Map all sub-shapes of the input solids
  TopExp::MapShapes(myInputShape, myInputsMap);

  // Collect all faces of the input shape requested for removal
  NCollection_List<TopoDS_Shape>           aFacesToRemove;
  NCollection_List<TopoDS_Shape>::Iterator aIt(myFacesToRemove);
  Message_ProgressScope                    aPSOuter(theRange, "Preparing features", 2);
  Message_ProgressScope                    aPS(aPSOuter.Next(),
                            "Preparing the faces to remove",
                            myFacesToRemove.Extent());
  for (; aIt.More(); aIt.Next(), aPS.Next())
  {
    const TopoDS_Shape& aS = aIt.Value();
    TopExp_Explorer     anExpF(aS, TopAbs_FACE);
    for (; anExpF.More(); anExpF.Next())
    {
      if (UserBreak(aPS))
      {
        return;
      }
      const TopoDS_Shape& aF = anExpF.Current();
      if (myInputsMap.Contains(aF))
      {
        aFacesToRemove.Append(aF);
      }
    }
  }

  if (aFacesToRemove.IsEmpty())
  {
    // Add error, that no features to remove have been found
    AddError(new BOPAlgo_AlertNoFacesToRemove());
    return;
  }

  // Build connexity blocks of the faces to remove
  TopoDS_Compound aCFToRemove;
  BRep_Builder().MakeCompound(aCFToRemove);
  for (aIt.Initialize(aFacesToRemove); aIt.More(); aIt.Next())
  {
    BRep_Builder().Add(aCFToRemove, aIt.Value());
  }

  // Fill the list of features with connexity blocks of faces
  BOPTools_AlgoTools::MakeConnexityBlocks(aCFToRemove, TopAbs_EDGE, TopAbs_FACE, myFeatures);
}

#include "BOPAlgo_RemoveFeatures_FillGap.pxx"

//=================================================================================================

//! Replaces the faces of the current shape by the faces of the input shape they have been made
//! from by the removal of the previous features
static TopoDS_Compound OriginalFaces(
  const NCollection_List<TopoDS_Shape>&                                            theFaces,
  const NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher>& theOrigins)
{
  TopoDS_Compound aCFaces;
  BRep_Builder().MakeCompound(aCFaces);
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aFence;
  NCollection_List<TopoDS_Shape>::Iterator               itLF(theFaces);
  for (; itLF.More(); itLF.Next())
  {
    const TopoDS_Shape* pOrigin = theOrigins.Seek(itLF.Value());
    const TopoDS_Shape& aF      = pOrigin ? *pOrigin : itLF.Value();
    if (aFence.Add(aF))
    {
      BRep_Builder().Add(aCFaces, aF);
    }
  }
  return aCFaces;
}

//=================================================================================================

//! Adds the faces made from the faces of the input feature by the removal of the previous
//! features (the faces themselves if they have not been modified)
static void CurrentFaces(const TopoDS_Shape&                                            theFeature,
                         const occ::handle<BRepTools_History>&                          theHistory,
                         NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>& theFaces)
{
  for (TopoDS_Iterator itF(theFeature); itF.More(); itF.Next())
  {
    const NCollection_List<TopoDS_Shape>& aLFIm = theHistory->Modified(itF.Value());
    if (aLFIm.IsEmpty())
    {
      theFaces.Add(itF.Value());
      continue;
    }
    for (NCollection_List<TopoDS_Shape>::Iterator itLFIm(aLFIm); itLFIm.More(); itLFIm.Next())
    {
      theFaces.Add(itLFIm.Value());
    }
  }
}

//=================================================================================================

//! Merges the faces and edges split by the removal of a feature, so that the next features find
//! whole adjacent faces. The sub-shapes of the shape the feature has been removed from are kept,
//! and so are the edges of the faces of the features still to be removed: merged with other
//! faces, they would take them away with the feature.
static TopoDS_Shape UnifyRemoval(
  const TopoDS_Shape&                                                  theShape,
  const TopoDS_Shape&                                                  thePrevShape,
  const NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>& theKeptFaces,
  const occ::handle<BRepTools_History>&                                theHistory)
{
  ShapeUpgrade_UnifySameDomain aSDTool;
  aSDTool.Initialize(theShape, true, true);
  aSDTool.AllowInternalEdges(false);
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aMS, aMPrev;
  TopExp::MapShapes(theShape, aMS);
  TopExp::MapShapes(thePrevShape, aMPrev);
  for (int i = 1; i <= aMPrev.Extent(); ++i)
  {
    if (aMS.Contains(aMPrev(i)))
    {
      aSDTool.KeepShape(aMPrev(i));
    }
  }
  for (int i = 1; i <= theKeptFaces.Extent(); ++i)
  {
    for (TopExp_Explorer anExpE(theKeptFaces(i), TopAbs_EDGE); anExpE.More(); anExpE.Next())
    {
      if (aMS.Contains(anExpE.Current()))
      {
        aSDTool.KeepShape(anExpE.Current());
      }
    }
  }
  try
  {
    OCC_CATCH_SIGNALS
    aSDTool.Build();
  }
  catch (Standard_Failure const&)
  {
    // The faces stay split, which only costs the next features some time
    return theShape;
  }
  theHistory->Merge(aSDTool.History());
  return aSDTool.Shape();
}

//=======================================================================
// function: RemoveFeatures
// purpose: Remove features by filling the gaps by extension of the
//          adjacent faces
//=======================================================================
void BOPAlgo_RemoveFeatures::RemoveFeatures(const Message_ProgressRange& theRange)
{
  // For each feature:
  // - Find the faces adjacent to the feature;
  // - Extend the adjacent faces;
  // - Trim the extended faces to fill the gap;
  // - Rebuild the solids with reconstructed adjacent faces
  //   and avoiding the feature faces.
  //
  // The gap of each feature is filled in the shape the previous features have been removed
  // from. Features often share adjacent faces (holes in one wall), and the faces rebuilt from
  // the input shape would bring back the gaps of the features removed before; the solids
  // would then have to be rebuilt from all of their faces, overlapping the rebuilt face with
  // its image from the previous removal, which is slow and does not always close the gaps
  // the same way.

  // Even if the history is not requested, it is necessary to track the faces
  // modification after each feature removal to find the faces of the next
  // features and to report the faces of the input shape.
  if (myHistory.IsNull())
  {
    myHistory = new BRepTools_History();
  }

  // Edge-Face connection map of the current shape to find faces adjacent to the feature,
  // and Face-Solid connection map to find the solids participating in its removal.
  // Rebuilt after each removed feature.
  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    anEFMap, anFSMap;
  // The faces of the input shape the faces of the current shape have been made from
  NCollection_DataMap<TopoDS_Shape, TopoDS_Shape, TopTools_ShapeMapHasher> anOrigins;
  bool isMapsValid = false;

  // Remove the features one by one.
  // It will allow removing the features even if there were
  // some problems with removal of the previous features.
  const int             aNbF = myFeatures.Extent();
  Message_ProgressScope aPSLoop(theRange, "Removing features one by one", aNbF);
  NCollection_List<TopoDS_Shape>::Iterator itF(myFeatures);
  for (; itF.More(); itF.Next())
  {
    if (UserBreak(aPSLoop))
    {
      return;
    }
    const TopoDS_Shape&   anInputFeature = itF.Value();
    Message_ProgressScope aPSFeature(aPSLoop.Next(), nullptr, 6);

    if (!isMapsValid)
    {
      anEFMap.Clear();
      anFSMap.Clear();
      TopExp::MapShapesAndAncestors(myShape, TopAbs_EDGE, TopAbs_FACE, anEFMap);
      TopExp::MapShapesAndAncestors(myShape, TopAbs_FACE, TopAbs_SOLID, anFSMap);
      anOrigins.Clear();
      for (int i = 1; i <= myInputsMap.Extent(); ++i)
      {
        if (myInputsMap(i).ShapeType() != TopAbs_FACE)
        {
          continue;
        }
        NCollection_List<TopoDS_Shape>::Iterator itLIm(myHistory->Modified(myInputsMap(i)));
        for (; itLIm.More(); itLIm.Next())
        {
          anOrigins.Bind(itLIm.Value(), myInputsMap(i));
        }
      }
      isMapsValid = true;
    }

    // The faces of the feature in the current shape
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aFeatureFaces;
    CurrentFaces(anInputFeature, myHistory, aFeatureFaces);
    TopoDS_Compound aFeature;
    BRep_Builder().MakeCompound(aFeature);
    bool bHasFaces = false;
    for (int i = 1; i <= aFeatureFaces.Extent(); ++i)
    {
      if (anFSMap.Contains(aFeatureFaces(i)))
      {
        BRep_Builder().Add(aFeature, aFeatureFaces(i));
        bHasFaces = true;
      }
    }
    if (!bHasFaces)
    {
      // Nothing is left of the feature
      continue;
    }

    // After the removal, the faces split by it are merged for the next features (the last
    // one leaves it to the simplification of the result), keeping the faces of the features
    // still to be removed
    const TopoDS_Shape aPrevShape = myShape;
    const auto         aOnRemoved = [&]() {
      isMapsValid = false;
      NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aRemainingFaces;
      NCollection_List<TopoDS_Shape>::Iterator                      itR(itF);
      for (itR.Next(); itR.More(); itR.Next())
      {
        CurrentFaces(itR.Value(), myHistory, aRemainingFaces);
      }
      if (!aRemainingFaces.IsEmpty())
      {
        myShape = UnifyRemoval(myShape, aPrevShape, aRemainingFaces, myHistory);
      }
    };
    FillGap aFG;
    aFG.SetFeature(aFeature);
    aFG.SetEFConnectionMap(anEFMap);
    aFG.SetFSConnectionMap(anFSMap);
    aFG.SetRunParallel(myRunParallel);
    aFG.SetFuzzyValue(myFuzzyValue);
    aFG.SetOriginalFaces(anOrigins);
    aFG.SetToSkipGrazing(mySkipGrazingFeatures);
    aFG.SetRange(aPSFeature.Next());
    aFG.Perform();
    if (UserBreak(aPSFeature))
    {
      return;
    }

    // A feature with grazing adjacent faces which could not be removed: skipped, or the
    // patch covering the gap next to them has not closed the solids
    const auto aReportGrazing = [&]() {
      AddWarning(new BOPAlgo_AlertGrazingAdjacentFaces(OriginalFaces(aFG.GrazingFaces(), anOrigins)));
      AddWarning(new BOPAlgo_AlertUnableToRemoveTheFeature(anInputFeature));
    };
    if (aFG.IsGrazing())
    {
      aReportGrazing();
      continue;
    }

    // Perform removal of the single feature
    if (RemoveFeature(aFG.Solids(),
                      aFG.FeatureFacesMap(),
                      aFG.HasAdjacentFaces(),
                      aFG.Faces(),
                      aFG.AnchoredFaces(),
                      aFG.History(),
                      HasHistory(),
                      aPSFeature.Next()))
    {
      aOnRemoved();
      if (aFG.IsPatched())
      {
        // Checked in the final result, the simplification may merge its edges
        myPatches.Append(aFG.Patch());
        myPatchSizes.Append(aFG.FeatureSize());
      }
      continue;
    }
    if (UserBreak(aPSFeature))
    {
      return;
    }
    // The extension of the grazing faces, which the patch has replaced, is not tried
    if (aFG.IsPatched())
    {
      aReportGrazing();
      continue;
    }

    // The adjacent faces extended by the size of the whole feature have not closed the
    // solid. The gap next to an adjacent face is as wide as the part of the feature the
    // face borders, which is much smaller than the whole feature for long features (ribs,
    // grooves) and features made of parts of different sizes; the faces extended that far
    // interfere with each other away from the gap. Retry the feature with every face
    // extended by the size of the feature faces it borders.
    FillGap aFGLocal;
    aFGLocal.SetFeature(aFeature);
    aFGLocal.SetEFConnectionMap(anEFMap);
    aFGLocal.SetFSConnectionMap(anFSMap);
    aFGLocal.SetRunParallel(myRunParallel);
    aFGLocal.SetFuzzyValue(myFuzzyValue);
    aFGLocal.SetOriginalFaces(anOrigins);
    aFGLocal.SetLocalExtension(true);
    aFGLocal.SetRange(aPSFeature.Next());
    aFGLocal.Perform();
    if (UserBreak(aPSFeature))
    {
      return;
    }
    if (RemoveFeature(aFGLocal.Solids(),
                      aFGLocal.FeatureFacesMap(),
                      aFGLocal.HasAdjacentFaces(),
                      aFGLocal.Faces(),
                      aFGLocal.AnchoredFaces(),
                      aFGLocal.History(),
                      HasHistory(),
                      aPSFeature.Next()))
    {
      aOnRemoved();
      continue;
    }
    if (UserBreak(aPSFeature))
    {
      return;
    }

    // The extensions keep only the tangent of a B-spline face, which runs off the surface it
    // was cut from: a wall split by the feature into pieces of one surface is not put back
    // together. Retry with the spans of all B-spline faces continued and such pieces rebuilt
    // as one face; tried only after the other attempts have failed, since the continuation
    // breaks some features they remove.
    FillGap aFGSpline;
    aFGSpline.SetFeature(aFeature);
    aFGSpline.SetEFConnectionMap(anEFMap);
    aFGSpline.SetFSConnectionMap(anFSMap);
    aFGSpline.SetRunParallel(myRunParallel);
    aFGSpline.SetFuzzyValue(myFuzzyValue);
    aFGSpline.SetOriginalFaces(anOrigins);
    aFGSpline.SetToSkipGrazing(mySkipGrazingFeatures);
    aFGSpline.SetSplineContinuation(true);
    aFGSpline.SetRange(aPSFeature.Next());
    aFGSpline.Perform();
    if (UserBreak(aPSFeature))
    {
      return;
    }
    if (!aFGSpline.IsGrazing()
        && RemoveFeature(aFGSpline.Solids(),
                         aFGSpline.FeatureFacesMap(),
                         aFGSpline.HasAdjacentFaces(),
                         aFGSpline.Faces(),
                         aFGSpline.AnchoredFaces(),
                         aFGSpline.History(),
                         HasHistory(),
                         aPSFeature.Next()))
    {
      aOnRemoved();
      continue;
    }
    if (UserBreak(aPSFeature))
    {
      return;
    }

    // Report the failure with the faces the first reconstruction could not rebuild
    AddWarning(new BOPAlgo_AlertUnableToRemoveTheFeature(anInputFeature));
    const TopoDS_Compound aNotRebuilt = OriginalFaces(aFG.FacesNotRebuilt(), anOrigins);
    for (TopoDS_Iterator itNR(aNotRebuilt); itNR.More(); itNR.Next())
    {
      AddWarning(new BOPAlgo_AlertUnableToRebuildAdjacentFace(itNR.Value()));
    }
  }
}

//=================================================================================================

//! Greatest tolerance of the edges and vertices of a patch, relative to the size of the feature
static constexpr double THE_PATCH_MAX_TOLERANCE = 1.e-3;

void BOPAlgo_RemoveFeatures::CheckPatchTolerance(const TopoDS_Shape& thePatch,
                                                 const double        theFeatureSize)
{
  if (!HasHistory())
  {
    // The pieces of the patch cannot be found in the simplified result
    return;
  }
  // The pieces of the patch in the result: the corners where the patch meets two faces
  // around it at small angles to each other may get far apart, and the vertices there take
  // up the distance in their tolerance
  const double   aMaxTol = THE_PATCH_MAX_TOLERANCE * theFeatureSize;
  TopoDS_Compound aCExceeding;
  BRep_Builder().MakeCompound(aCExceeding);
  bool                                                   isExceeding = false;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aMChecked;
  // The history follows the patch itself, its pieces are its modifications
  NCollection_List<TopoDS_Shape>::Iterator itF(myHistory->Modified(thePatch));
  for (; itF.More(); itF.Next())
  {
    for (TopExp_Explorer anExp(itF.Value(), TopAbs_EDGE); anExp.More(); anExp.Next())
    {
      if (aMChecked.Add(anExp.Current())
          && BRep_Tool::Tolerance(TopoDS::Edge(anExp.Current())) > aMaxTol)
      {
        BRep_Builder().Add(aCExceeding, anExp.Current());
        isExceeding = true;
      }
    }
    for (TopExp_Explorer anExp(itF.Value(), TopAbs_VERTEX); anExp.More(); anExp.Next())
    {
      if (aMChecked.Add(anExp.Current())
          && BRep_Tool::Tolerance(TopoDS::Vertex(anExp.Current())) > aMaxTol)
      {
        BRep_Builder().Add(aCExceeding, anExp.Current());
        isExceeding = true;
      }
    }
  }
  if (isExceeding)
  {
    AddWarning(new BOPAlgo_AlertPatchTolerance(aCExceeding));
  }
}

//=================================================================================================

//! An image of an adjacent face narrower than this fraction of the size of the feature may
//! vanish in the reconstruction of the solid: the intersection collapses such slivers
static constexpr double THE_SLIVER_WIDTH = 1.e-5;

//! Checks whether every face of the shape is narrower than the given width (its area is
//! less than the width times the size of its box)
static bool IsSliver(const TopoDS_Shape& theShape, const double theMaxWidth)
{
  for (TopExp_Explorer anExp(theShape, TopAbs_FACE); anExp.More(); anExp.Next())
  {
    Bnd_Box aBox;
    BRepBndLib::AddOptimal(anExp.Current(), aBox, false, false);
    GProp_GProps aProps;
    BRepGProp::SurfaceProperties(anExp.Current(), aProps);
    if (!aBox.IsVoid() && aProps.Mass() >= theMaxWidth * sqrt(aBox.SquareExtent()))
    {
      return false;
    }
  }
  return true;
}

//=================================================================================================

bool BOPAlgo_RemoveFeatures::RemoveFeature(
  const NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>& theSolids,
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>&        theFeatureFacesMap,
  const bool                                                           theHasAdjacentFaces,
  const NCollection_IndexedDataMap<TopoDS_Shape,
                                   NCollection_List<TopoDS_Shape>,
                                   TopTools_ShapeMapHasher>&           theAdjFaces,
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>&        theAnchoredFaces,
  const occ::handle<BRepTools_History>&                                theAdjFacesHistory,
  const bool                                                           theSolidsHistoryNeeded,
  const Message_ProgressRange&                                         theRange)
{
  bool      bFuseShapes = true;
  const int aNbAF       = theAdjFaces.Extent();
  if (aNbAF == 0)
  {
    if (theHasAdjacentFaces)
    {
      // The adjacent faces have been found for the feature,
      // but something went wrong during their rebuilding.
      return false;
    }

    // No adjacent faces have been found for the feature.
    // Most likely, the intention is to remove some internal part of the shape.
    // We just have to rebuild the shape without the feature, no fuse is necessary.
    bFuseShapes = false;
  }

  Message_ProgressScope aPS(theRange, nullptr, 100);
  // Rebuild the shape using MakerVolume algorithm avoiding the faces of the
  // feature and replacing the adjacent faces with their images

  BRep_Builder aBB;

  // From the faces of input shape build map of faces with which the result will be validated:
  // - all non-internal faces of the input shape (except for adjacent and feature faces)
  //   must have some trace in the result solids;
  // - all adjacent faces (if any) must have some trace in the result solids.
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aFacesToBeKept;
  // Build also the map of faces to compare orientation of the original and result faces
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aFacesToCheckOri;

  // Collect list of internal entities of the input shape to be avoided in result
  // and to make them removed in the history.
  NCollection_List<TopoDS_Shape> anInternalShapes;

  // In the feature removal will participate only the solids connected to the feature
  // or the faces adjacent to the feature.

  // Solids to be rebuilt  after the feature removal
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aSolidsToRebuild;
  // Find faces shared between solids
  NCollection_List<TopoDS_Shape> aSharedFaces;
  // Solids to be avoided in the feature removal and added into result directly
  NCollection_List<TopoDS_Shape> anUnTouchedSolids;

  // Prepare to the feature removal - fill all necessary containers
  GetOriginalFaces(myShape,
                   theSolids,
                   theFeatureFacesMap,
                   theAdjFaces,
                   myHistory,
                   aFacesToBeKept,
                   anInternalShapes,
                   aFacesToCheckOri,
                   aSolidsToRebuild,
                   aSharedFaces,
                   anUnTouchedSolids);

  aPS.Next(3);

  // Usually the reconstructed adjacent faces fit into the gap left by the feature exactly, and
  // the solid can be rebuilt by replacing them in its shell. That only needs the faces around
  // them to be intersected with them, while the general reconstruction below intersects all
  // faces of the solid and builds the solid from scratch, which takes most of the time of a
  // small feature on a big model. Whenever the faces do not fit exactly, the general
  // reconstruction takes over.
  if (bFuseShapes && aSolidsToRebuild.Extent() == 1 && aSharedFaces.IsEmpty()
      && anInternalShapes.IsEmpty())
  {
    BOPAlgo_Builder    aGFLocal;
    const TopoDS_Shape aSolid = RebuildSolidLocally(aSolidsToRebuild(1),
                                                    theFeatureFacesMap,
                                                    theAdjFaces,
                                                    theAdjFacesHistory,
                                                    myFuzzyValue,
                                                    myRunParallel,
                                                    aGFLocal,
                                                    aPS.Next(20));
    if (!aSolid.IsNull())
    {
      myHistory->Merge(theAdjFacesHistory);
      myHistory->Merge(aGFLocal.History());
      if (theSolidsHistoryNeeded)
      {
        BRepTools_History aSolidsHistory;
        aSolidsHistory.AddModified(aSolidsToRebuild(1), aSolid);
        myHistory->Merge(aSolidsHistory);
      }

      TopoDS_Compound aCRes;
      aBB.MakeCompound(aCRes);
      aBB.Add(aCRes, aSolid);
      NCollection_List<TopoDS_Shape>::Iterator itLS(anUnTouchedSolids);
      for (; itLS.More(); itLS.Next())
      {
        aBB.Add(aCRes, itLS.Value());
      }
      myShape = aCRes;
      return true;
    }
    if (UserBreak(aPS))
    {
      return false;
    }
  }

  // To avoid excessive intersection of the faces collect the faces
  // of the input shape into a compound
  TopoDS_Compound anOrigF;
  aBB.MakeCompound(anOrigF);
  int aNbFK = aFacesToBeKept.Extent();
  for (int i = 1; i <= aNbFK; ++i)
  {
    aBB.Add(anOrigF, aFacesToBeKept(i));
  }

  // Tool for solids reconstruction
  BOPAlgo_MakerVolume aMV;
  aMV.SetRunParallel(myRunParallel);
  aMV.SetFuzzyValue(myFuzzyValue);
  aMV.SetAvoidInternalShapes(true);
  aMV.SetIntersect(bFuseShapes);
  aMV.SetNonDestructive(true);
  // Add faces of the input shape
  aMV.AddArgument(anOrigF);

  // Add reconstructed adjacent faces. The faces rebuilt as a group of coincident
  // faces share their images, which are added once.
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> anAdjImagesAdded;
  for (int i = 1; i <= aNbAF; ++i)
  {
    if (HasHistory())
    {
      // Look for internal edges in the original adjacent faces
      const TopoDS_Shape& aFOr = theAdjFaces.FindKey(i);
      FindInternals(aFOr, anInternalShapes);
    }

    NCollection_List<TopoDS_Shape>           aLFA;
    NCollection_List<TopoDS_Shape>::Iterator itLFAIm(theAdjFaces(i));
    for (; itLFAIm.More(); itLFAIm.Next())
    {
      if (anAdjImagesAdded.Add(itLFAIm.Value()))
      {
        aLFA.Append(itLFAIm.Value());
      }
    }
    if (aLFA.IsEmpty())
    {
      // The images have already been added with another face of the group
      continue;
    }

    if (aLFA.Extent() == 1)
    {
      const TopoDS_Shape& aFA = aLFA.First();
      aMV.AddArgument(aFA);
      aFacesToBeKept.Add(aFA);
    }
    else
    {
      // To avoid intersection among the images, collect them into compound
      TopoDS_Compound anAdjF;
      aBB.MakeCompound(anAdjF);
      NCollection_List<TopoDS_Shape>::Iterator itLFA(aLFA);
      for (; itLFA.More(); itLFA.Next())
      {
        aBB.Add(anAdjF, itLFA.Value());
      }

      aMV.AddArgument(anAdjF);
      aFacesToBeKept.Add(anAdjF);
    }
  }

  // Build solids
  aMV.Perform(aPS.Next(70));
  if (aMV.HasErrors())
  {
    return false;
  }

  // Result of MV operation
  const TopoDS_Shape& aSolids = aMV.Shape();
  TopExp_Explorer     anExpS(aSolids, TopAbs_SOLID);
  if (!anExpS.More())
  {
    // No solids have been built
    return false;
  }

  // Now, it is necessary to:
  // - Validate the result by checking faces of the map <aFacesToBeKept>
  //   to have some parts kept in the resulting solids;
  // - Remove the solids possibly filling the holes in the original shapes;
  // - Update history with the history of MakerVolume algorithm.

  // Splits of adjacent faces from previous runs
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> anAdjFacesSplits;
  for (int i = 1; i <= aNbAF; ++i)
  {
    if (!aPS.More())
    {
      return false;
    }
    const TopoDS_Shape&                   aF    = theAdjFaces.FindKey(i);
    const NCollection_List<TopoDS_Shape>& aLFIm = myHistory->Modified(aF);
    if (aLFIm.IsEmpty())
    {
      anAdjFacesSplits.Add(aF);
    }
    else
    {
      NCollection_List<TopoDS_Shape>::Iterator itLFIm(aLFIm);
      for (; itLFIm.More(); itLFIm.Next())
      {
        anAdjFacesSplits.Add(itLFIm.Value());
      }
    }
  }

  // Validate the result
  Bnd_Box aFeatureBox;
  for (const TopoDS_Shape& aFF : theFeatureFacesMap)
  {
    BRepBndLib::Add(aFF, aFeatureBox);
  }
  const double aSliverWidth =
    aFeatureBox.IsVoid() ? 0. : THE_SLIVER_WIDTH * sqrt(aFeatureBox.SquareExtent());
  bool bValid = true;
  aNbFK       = aFacesToBeKept.Extent();
  for (int i = 1; i <= aNbFK && bValid; ++i)
  {

    const TopoDS_Shape& aS = aFacesToBeKept(i);
    if (anAdjFacesSplits.Contains(aS))
    {
      continue;
    }
    if (!aPS.More())
    {
      return false;
    }
    TopExp_Explorer anExpF(aS, TopAbs_FACE);
    for (; anExpF.More(); anExpF.Next())
    {
      const TopoDS_Shape& aF = anExpF.Current();
      if (!aMV.IsDeleted(aF))
      {
        break;
      }
    }
    bValid = anExpF.More()
             || (anAdjImagesAdded.Contains(TopExp_Explorer(aS, TopAbs_FACE).Current())
                 && IsSliver(aS, aSliverWidth));
  }

  if (!bValid)
  {
    return false;
  }

  // It is possible that the result of MakerVolume operation contains
  // some extra solids. Get only correct solids.
  NCollection_List<TopoDS_Shape> aLSRes;
  // Remember the removed shapes
  NCollection_List<TopoDS_Shape> aRemovedShapes;
  GetValidSolids(aMV,
                 aFacesToCheckOri,
                 theAdjFaces,
                 theAnchoredFaces,
                 aSharedFaces,
                 anOrigF,
                 theSolids.Extent(),
                 aLSRes,
                 aRemovedShapes);

  if (aLSRes.Extent() != theSolids.Extent())
  {
    return false;
  }
  aPS.Next(3);
  // Remove internal wires from the faces, possibly appeared after intersection
  RemoveInternalWires(aLSRes, &anInternalShapes);

  // Update history with:
  // History of adjacent faces reconstruction
  myHistory->Merge(theAdjFacesHistory);
  // History of intersection
  myHistory->Merge(aMV.History());

  if (HasHistory())
  {
    // Map the result to check if the shape is removed
    NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aMSRes;
    NCollection_List<TopoDS_Shape>::Iterator                      itLS(aLSRes);
    for (; itLS.More(); itLS.Next())
    {
      TopExp::MapShapes(itLS.Value(), aMSRes);
    }

    // Remove internal shapes and extra faces
    BRepTools_History aRemHist;
    anInternalShapes.Append(aRemovedShapes);
    MakeRemoved(anInternalShapes, aRemHist, aMSRes);
    myHistory->Merge(aRemHist);
  }
  aPS.Next(3);
  // Fill the history for the solids
  if (theSolidsHistoryNeeded)
  {
    BRepTools_History aSolidsHistory;
    FillSolidsHistory(aSolidsToRebuild, aLSRes, theAdjFaces, aMV, aSolidsHistory);
    myHistory->Merge(aSolidsHistory);
  }

  TopoDS_Compound aCRes;
  aBB.MakeCompound(aCRes);
  // Add reconstructed solids
  NCollection_List<TopoDS_Shape>::Iterator itLS(aLSRes);
  for (; itLS.More(); itLS.Next())
  {
    aBB.Add(aCRes, itLS.Value());
  }

  // Add unmodified solids
  itLS.Initialize(anUnTouchedSolids);
  for (; itLS.More(); itLS.Next())
  {
    aBB.Add(aCRes, itLS.Value());
  }

  // Save the result
  myShape = aCRes;
  return true;
}

//=======================================================================
// function: UpdateHistory
// purpose: Update history with the removed features
//=======================================================================
void BOPAlgo_RemoveFeatures::UpdateHistory(const Message_ProgressRange& theRange)
{
  if (!HasHistory())
  {
    return;
  }

  // Map the result
  myMapShape.Clear();
  TopExp::MapShapes(myShape, myMapShape);

  // Update the history
  BRepTools_History aHistory;

  const int             aNbS = myInputsMap.Extent();
  Message_ProgressScope aPS(theRange, "Updating history", aNbS);
  for (int i = 1; i <= aNbS; ++i, aPS.Next())
  {
    const TopoDS_Shape& aS = myInputsMap(i);
    if (!BRepTools_History::IsSupportedType(aS))
    {
      continue;
    }

    if (myHistory->IsRemoved(aS))
    {
      continue;
    }

    if (UserBreak(aPS))
    {
      return;
    }
    // Check if the shape has any trace in the result
    const NCollection_List<TopoDS_Shape>& aLSIm = myHistory->Modified(aS);
    if (aLSIm.IsEmpty())
    {
      if (!myMapShape.Contains(aS))
      {
        aHistory.Remove(aS);
      }
    }

    NCollection_List<TopoDS_Shape>::Iterator itLSIm(aLSIm);
    for (; itLSIm.More(); itLSIm.Next())
    {
      if (!myMapShape.Contains(itLSIm.Value()))
      {
        aHistory.Remove(itLSIm.Value());
      }
    }
  }

  myHistory->Merge(aHistory);
}

//=======================================================================
// function: SimplifyResult
// purpose: Simplifies the result by removing extra edges and vertices
//          created during operation
//=======================================================================
void BOPAlgo_RemoveFeatures::SimplifyResult(const Message_ProgressRange& theRange)
{
  if (myShape.IsSame(myInputShape))
  {
    return;
  }
  Message_ProgressScope        aPSOuter(theRange, "Simplifyingthe result", 2);
  ShapeUpgrade_UnifySameDomain aSDTool;
  aSDTool.Initialize(myShape, true, true);
  // Do not allow producing internal edges
  aSDTool.AllowInternalEdges(false);
  // Avoid removal of the input edges and vertices
  if (myMapShape.IsEmpty())
  {
    TopExp::MapShapes(myShape, myMapShape);
  }

  const int             aNbS = myInputsMap.Extent();
  Message_ProgressScope aPS(aPSOuter.Next(), nullptr, aNbS);
  for (int i = 1; i <= aNbS; ++i, aPS.Next())
  {
    if (UserBreak(aPS))
    {
      return;
    }
    if (myMapShape.Contains(myInputsMap(i)))
    {
      aSDTool.KeepShape(myInputsMap(i));
    }
  }

  // Perform unification. The result is valid without it, and is kept as it is if the
  // unification fails.
  try
  {
    OCC_CATCH_SIGNALS
    aSDTool.Build();
  }
  catch (Standard_Failure const&)
  {
    return;
  }
  aPSOuter.Next();
  myShape = aSDTool.Shape();
  if (HasHistory())
  {
    myHistory->Merge(aSDTool.History());
  }
}

//=======================================================================
// function: PostTreat
// purpose: Restore the type of the initial shape
//=======================================================================
void BOPAlgo_RemoveFeatures::PostTreat()
{
  const TopAbs_ShapeEnum anInputType = myInputShape.ShapeType();
  const TopAbs_ShapeEnum aResType    = myShape.ShapeType();
  if (aResType == anInputType)
  {
    return;
  }

  TopExp_Explorer anExpS(myShape, TopAbs_SOLID);

  if (anInputType == TopAbs_SOLID)
  {
    myShape = anExpS.Current();
    return;
  }

  TopoDS_Shape aRes;
  if (anInputType == TopAbs_COMPOUND)
  {
    BRep_Builder().MakeCompound(TopoDS::Compound(aRes));
  }
  else
  {
    BRep_Builder().MakeCompSolid(TopoDS::CompSolid(aRes));
  }

  for (; anExpS.More(); anExpS.Next())
  {
    BRep_Builder().Add(aRes, anExpS.Current());
  }

  myShape = aRes;
}

//=======================================================================
// static methods definition
//=======================================================================

//=======================================================================
// function: MakeRemoved
// purpose: Makes the shapes in the list removed in the history.
//          Keeps the shapes contained in the map.
//=======================================================================
void MakeRemoved(const NCollection_List<TopoDS_Shape>&                                theShapes,
                 BRepTools_History&                                                   theHistory,
                 const NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>& theKeepShapes)
{
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aShapesMap;
  NCollection_List<TopoDS_Shape>::Iterator                      it(theShapes);
  for (; it.More(); it.Next())
  {
    TopExp::MapShapes(it.Value(), aShapesMap);
  }

  const int aNbS = aShapesMap.Extent();
  for (int i = 1; i <= aNbS; ++i)
  {
    const TopoDS_Shape& aS = aShapesMap(i);
    if (!theKeepShapes.Contains(aS) && BRepTools_History::IsSupportedType(aS))
    {
      theHistory.Remove(aS);
    }
  }
}

//=======================================================================
// function: FindInternals
// purpose: Looks for internal shapes inside the face or solid
//=======================================================================
void FindInternals(const TopoDS_Shape& theS, NCollection_List<TopoDS_Shape>& theLInt)
{
  TopoDS_Iterator itS(theS);
  for (; itS.More(); itS.Next())
  {
    const TopoDS_Shape& aSS = itS.Value();
    if (aSS.Orientation() == TopAbs_INTERNAL)
    {
      theLInt.Append(aSS);
    }
    else
    {
      TopoDS_Iterator itSS(aSS);
      for (; itSS.More(); itSS.Next())
      {
        if (itSS.Value().Orientation() == TopAbs_INTERNAL)
        {
          theLInt.Append(aSS);
          break;
        }
      }
    }
  }
}

//=======================================================================
// function: RemoveInternalWires
// purpose: Removes internal wires from the faces
//=======================================================================
void RemoveInternalWires(const NCollection_List<TopoDS_Shape>& theShapes,
                         NCollection_List<TopoDS_Shape>*       theRemoved)
{
  NCollection_List<TopoDS_Shape>::Iterator itLS(theShapes);
  for (; itLS.More(); itLS.Next())
  {
    const TopoDS_Shape& aShape = itLS.Value();
    TopExp_Explorer     anExpF(aShape, TopAbs_FACE);
    for (; anExpF.More(); anExpF.Next())
    {
      TopoDS_Face&                   aF = *(TopoDS_Face*)&anExpF.Current();
      NCollection_List<TopoDS_Shape> aLWToRemove;
      FindInternals(aF, aLWToRemove);
      if (aLWToRemove.Extent())
      {
        aF.Free(true);
        NCollection_List<TopoDS_Shape>::Iterator itR(aLWToRemove);
        for (; itR.More(); itR.Next())
        {
          if (theRemoved)
          {
            theRemoved->Append(itR.Value());
          }
          BRep_Builder().Remove(aF, itR.Value());
        }
        aF.Free(false);
      }
    }
  }
}

//=======================================================================
// function: GetOriginalFaces
// purpose: Get original faces from my face
//=======================================================================
void GetOriginalFaces(
  const TopoDS_Shape&                                                  theShape,
  const NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>& theSolids,
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>&        theFeatureFacesMap,
  const NCollection_IndexedDataMap<TopoDS_Shape,
                                   NCollection_List<TopoDS_Shape>,
                                   TopTools_ShapeMapHasher>&           theAdjFaces,
  const occ::handle<BRepTools_History>&                                theHistory,
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>&       theFacesToBeKept,
  NCollection_List<TopoDS_Shape>&                                      theInternalShapes,
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>&              theFacesToCheckOri,
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>&       theSolidsToRebuild,
  NCollection_List<TopoDS_Shape>&                                      theSharedFaces,
  NCollection_List<TopoDS_Shape>&                                      theUnTouchedSolids)
{
  // Use only solids which has to be reconstructed by the feature.
  // All other solids should be avoided in the feature removal and added
  // into result directly.

  // Update solids after removal of previous features
  const int aNbSols = theSolids.Extent();
  for (int i = 1; i <= aNbSols; ++i)
  {
    const TopoDS_Shape&                   aSol  = theSolids(i);
    const NCollection_List<TopoDS_Shape>& aLFIm = theHistory->Modified(aSol);
    if (aLFIm.IsEmpty())
    {
      theSolidsToRebuild.Add(aSol);
    }
    else
    {
      theSolidsToRebuild.Add(aLFIm.First());
    }
  }

  // Splits of the feature faces
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>           aFeatureFacesSplits;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>::Iterator itM(theFeatureFacesMap);
  for (; itM.More(); itM.Next())
  {
    const TopoDS_Shape&                   aF    = itM.Value();
    const NCollection_List<TopoDS_Shape>& aLFIm = theHistory->Modified(aF);
    if (aLFIm.IsEmpty())
    {
      aFeatureFacesSplits.Add(aF);
    }
    else
    {
      NCollection_List<TopoDS_Shape>::Iterator itLFIm(aLFIm);
      for (; itLFIm.More(); itLFIm.Next())
      {
        aFeatureFacesSplits.Add(itLFIm.Value());
      }
    }
  }

  TopExp_Explorer anExpS(theShape, TopAbs_SOLID);
  for (; anExpS.More(); anExpS.Next())
  {
    const TopoDS_Shape& aSol = anExpS.Current();

    // Check if the solid has to be reconstructed
    if (!theSolidsToRebuild.Contains(aSol))
    {
      // untouched solid
      theUnTouchedSolids.Append(aSol);
      continue;
    }

    TopoDS_Iterator itSh(aSol);
    for (; itSh.More(); itSh.Next())
    {
      const TopoDS_Shape& aSh = itSh.Value();
      if (aSh.ShapeType() != TopAbs_SHELL)
      {
        theInternalShapes.Append(aSh);
        continue;
      }

      TopoDS_Iterator itF(aSh);
      for (; itF.More(); itF.Next())
      {
        const TopoDS_Shape& aF = itF.Value();
        // Avoid the feature faces
        if (aFeatureFacesSplits.Contains(aF))
        {
          continue;
        }

        // Avoid the adjacent faces
        if (theAdjFaces.Contains(aF))
        {
          continue;
        }

        if (aF.Orientation() != TopAbs_INTERNAL)
        {
          theFacesToBeKept.Add(aF);

          if (!theFacesToCheckOri.Add(aF))
          {
            theFacesToCheckOri.Remove(aF);
            theSharedFaces.Append(aF);
          }
        }
        else
        {
          theInternalShapes.Append(aSh);
        }
      }
    }
  }
}

//=======================================================================
// function: FindShape
// purpose: Find the shape in the other shape
//=======================================================================
void FindShape(const TopoDS_Shape& theSWhat, const TopoDS_Shape& theSWhere, TopoDS_Shape& theSFound)
{
  TopExp_Explorer anExp(theSWhere, theSWhat.ShapeType());
  for (; anExp.More(); anExp.Next())
  {
    const TopoDS_Shape& aS = anExp.Current();
    if (aS.IsSame(theSWhat))
    {
      theSFound = aS;
      break;
    }
  }
}

//=======================================================================
// function: GetValidSolids
// purpose: Checks the validity of the solids and keeps only valid ones
//=======================================================================
void GetValidSolids(
  BOPAlgo_MakerVolume&                                          theMV,
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theFacesToCheckOri,
  const NCollection_IndexedDataMap<TopoDS_Shape,
                                   NCollection_List<TopoDS_Shape>,
                                   TopTools_ShapeMapHasher>&    theAdjFaces,
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theAnchoredFaces,
  const NCollection_List<TopoDS_Shape>&                         aSharedFaces,
  const TopoDS_Shape&                                           theOrigFaces,
  const int                                                     theNbSol,
  NCollection_List<TopoDS_Shape>&                               theLSRes,
  NCollection_List<TopoDS_Shape>&                               theRemovedShapes)
{
  TopExp_Explorer anExpS(theMV.Shape(), TopAbs_SOLID);
  for (; anExpS.More(); anExpS.Next())
  {
    theLSRes.Append(anExpS.Current());
  }

  if (theLSRes.Extent() > theNbSol)
  {
    // Find Solids filling the holes in the initial shape
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aSolidsToAvoid;
    NCollection_IndexedDataMap<TopoDS_Shape,
                               NCollection_List<TopoDS_Shape>,
                               TopTools_ShapeMapHasher>
      aFSMap;
    TopExp::MapShapesAndAncestors(theMV.Shape(), TopAbs_FACE, TopAbs_SOLID, aFSMap);
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aValidSolids;
    FindExtraShapes(aFSMap, theFacesToCheckOri, theMV, aSolidsToAvoid, &aValidSolids);

    NCollection_List<TopoDS_Shape>::Iterator itLS(theLSRes);
    for (; itLS.More();)
    {
      if (aSolidsToAvoid.Contains(itLS.Value()))
      {
        theLSRes.Remove(itLS);
      }
      else
      {
        itLS.Next();
      }
    }

    // The solids containing no faces of the original shape are classified by the
    // reconstructed adjacent faces: a solid lying on the wrong side of any of them
    // is a void (e.g. the inside of a pocket whose mouth got closed by an extension
    // of the top face), unless a split validated by an edge of the original face
    // confirms the solid as material. Without this check such a void would be
    // merged into the result by AvoidExtraSharedFaces.
    if (theLSRes.Extent() > theNbSol)
    {
      NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aVoidSolids, aMaterialSolids;
      const int                                              aNbAF = theAdjFaces.Extent();
      for (int i = 1; i <= aNbAF; ++i)
      {
        NCollection_List<TopoDS_Shape>::Iterator itLFA(theAdjFaces(i));
        for (; itLFA.More(); itLFA.Next())
        {
          const TopoDS_Face&             aFA       = TopoDS::Face(itLFA.Value());
          const bool                     bAnchored = theAnchoredFaces.Contains(aFA);
          NCollection_List<TopoDS_Shape> aLFIm;
          TakeModified(aFA, theMV, aLFIm);
          NCollection_List<TopoDS_Shape>::Iterator itLFIm(aLFIm);
          for (; itLFIm.More(); itLFIm.Next())
          {
            const NCollection_List<TopoDS_Shape>* pLS = aFSMap.Seek(itLFIm.Value());
            if (!pLS)
            {
              continue;
            }
            NCollection_List<TopoDS_Shape>::Iterator itLS1(*pLS);
            for (; itLS1.More(); itLS1.Next())
            {
              const TopoDS_Shape& aSol = itLS1.Value();
              if (aValidSolids.Contains(aSol))
              {
                continue;
              }
              TopoDS_Face aFInSol;
              FindShape(itLFIm.Value(), aSol, aFInSol);
              if (aFInSol.IsNull())
              {
                continue;
              }
              // Compare with the reconstructed face, not with the original one: the
              // image usually lies outside the parametric bounds of the original face,
              // where the projection used by IsSplitToReverse is unreliable.
              if (BOPTools_AlgoTools::IsSplitToReverse(aFInSol, aFA, theMV.Context()))
              {
                aVoidSolids.Add(aSol);
              }
              else if (bAnchored)
              {
                aMaterialSolids.Add(aSol);
              }
            }
          }
        }
      }

      for (itLS.Initialize(theLSRes); itLS.More();)
      {
        if (aVoidSolids.Contains(itLS.Value()) && !aMaterialSolids.Contains(itLS.Value()))
        {
          theRemovedShapes.Append(itLS.Value());
          theLSRes.Remove(itLS);
        }
        else
        {
          itLS.Next();
        }
      }
    }
  }

  if (theLSRes.Extent() > theNbSol)
  {
    // Check if the splits of the adjacent faces split the solids
    AvoidExtraSharedFaces(theLSRes, aSharedFaces, theMV, theRemovedShapes);
  }

  if (theLSRes.Extent() > theNbSol)
  {
    // Remove solids containing only the adjacent faces
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> anOrigFacesRes;
    TopExp_Explorer                                        anExpF(theOrigFaces, TopAbs_FACE);
    for (; anExpF.More(); anExpF.Next())
    {
      TakeModified(anExpF.Current(), theMV, anOrigFacesRes);
    }

    NCollection_List<TopoDS_Shape>::Iterator itLS(theLSRes);
    for (; itLS.More();)
    {
      anExpF.Init(itLS.Value(), TopAbs_FACE);
      for (; anExpF.More(); anExpF.Next())
      {
        if (anOrigFacesRes.Contains(anExpF.Current()))
        {
          break;
        }
      }
      if (!anExpF.More())
      {
        theRemovedShapes.Append(itLS.Value());
        theLSRes.Remove(itLS);
      }
      else
      {
        itLS.Next();
      }
    }
  }
}

//=======================================================================
// function: FindExtraShape
// purpose: Find shapes possibly filling the holes in the original shape
//=======================================================================
void FindExtraShapes(
  const NCollection_IndexedDataMap<TopoDS_Shape,
                                   NCollection_List<TopoDS_Shape>,
                                   TopTools_ShapeMapHasher>&    theConnectionMap,
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theShapesToCheckOri,
  BOPAlgo_Builder&                                              theBuilder,
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>&       theShapesToAvoid,
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>*       theValidShapes,
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>* theInternalOK)
{
  occ::handle<IntTools_Context>                           aCtx = theBuilder.Context();
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>  aValidShapes;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>* pValidShapes =
    theValidShapes ? theValidShapes : &aValidShapes;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>::Iterator itM(theShapesToCheckOri);
  for (; itM.More(); itM.Next())
  {
    const TopoDS_Shape& aSToCheckOri = itM.Value();
    // Check modification of the shape during intersection
    NCollection_List<TopoDS_Shape> aLSIm;
    TakeModified(aSToCheckOri, theBuilder, aLSIm);

    NCollection_List<TopoDS_Shape>::Iterator itLSIm(aLSIm);
    for (; itLSIm.More(); itLSIm.Next())
    {
      const TopoDS_Shape& aSIm = itLSIm.Value();

      const NCollection_List<TopoDS_Shape>* pShapesToValidate = theConnectionMap.Seek(aSIm);
      if (!pShapesToValidate)
      {
        continue;
      }

      NCollection_List<TopoDS_Shape>::Iterator itSV(*pShapesToValidate);
      for (; itSV.More(); itSV.Next())
      {
        const TopoDS_Shape& aShapeToValidate = itSV.Value();
        if (pValidShapes->Contains(aShapeToValidate))
        {
          continue;
        }

        TopoDS_Face aSInShape;
        FindShape(aSIm, aShapeToValidate, aSInShape);

        // An edge that is INTERNAL in a split is not shared with any other split, so its
        // orientation says nothing about the split: the split spans both sides of the edge.
        // It is legitimate only when the faces on the other side of the edge are themselves
        // reconstructed or removed (see TrimFace), otherwise the split overlaps a kept face.
        bool bSameOri = (aSInShape.Orientation() == TopAbs_INTERNAL && theInternalOK)
                          ? theInternalOK->Contains(aSToCheckOri)
                          : !BOPTools_AlgoTools::IsSplitToReverse(aSInShape, aSToCheckOri, aCtx);

        if (bSameOri)
        {
          pValidShapes->Add(aShapeToValidate);
        }
        else
        {
          theShapesToAvoid.Add(aShapeToValidate);
        }
      }
    }
  }

  itM.Initialize(*pValidShapes);
  for (; itM.More(); itM.Next())
  {
    theShapesToAvoid.Remove(itM.Value());
  }
}

//=======================================================================
// function: AvoidExtraSharedFaces
// purpose: Looks for the extra faces splitting the solids
//=======================================================================
void AvoidExtraSharedFaces(NCollection_List<TopoDS_Shape>&       theLSolids,
                           const NCollection_List<TopoDS_Shape>& theLFSharedToAvoid,
                           BOPAlgo_Builder&                      theBuilder,
                           NCollection_List<TopoDS_Shape>&       theExtraFaces)
{
  // Get all splits of shared faces to avoid in the check
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aMFSharedSp;
  {
    NCollection_List<TopoDS_Shape>           aLFSharedSp;
    NCollection_List<TopoDS_Shape>::Iterator itLFS(theLFSharedToAvoid);
    for (; itLFS.More(); itLFS.Next())
    {
      TakeModified(itLFS.Value(), theBuilder, aMFSharedSp);
    }
  }

  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
                                           aFSMap;
  NCollection_List<TopoDS_Shape>::Iterator itLS(theLSolids);
  for (; itLS.More(); itLS.Next())
  {
    TopExp::MapShapesAndAncestors(itLS.Value(), TopAbs_FACE, TopAbs_SOLID, aFSMap);
  }

  NCollection_List<TopoDS_Shape> anExtraFaces;
  NCollection_List<TopoDS_Shape> aLFArguments;
  itLS.Initialize(theLSolids);
  for (; itLS.More(); itLS.Next())
  {
    const TopoDS_Shape& aSol = itLS.Value();
    TopExp_Explorer     anExpF(aSol, TopAbs_FACE);
    for (; anExpF.More(); anExpF.Next())
    {
      const TopoDS_Shape&                   aF    = anExpF.Current();
      const NCollection_List<TopoDS_Shape>& aLSol = aFSMap.FindFromKey(aF);
      if (aLSol.Extent() != 2 || aMFSharedSp.Contains(aF))
      {
        aLFArguments.Append(aF);
      }
      else
      {
        anExtraFaces.Append(aF);
      }
    }
  }

  if (anExtraFaces.IsEmpty())
  {
    return;
  }

  // Rebuild the solids avoiding the extra faces
  BOPAlgo_BuilderSolid aBS;
  aBS.SetAvoidInternalShapes(true);
  aBS.SetShapes(aLFArguments);
  aBS.Perform();
  if (aBS.HasErrors())
  {
    return;
  }

  theLSolids = aBS.Areas();
  theExtraFaces.Append(anExtraFaces);
}

//=======================================================================
// function: FillSolidsHistory
// purpose: Fills the history of solids modifications
//=======================================================================
void FillSolidsHistory(
  const NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher>& theSolIn,
  NCollection_List<TopoDS_Shape>&                                      theSolidsOut,
  const NCollection_IndexedDataMap<TopoDS_Shape,
                                   NCollection_List<TopoDS_Shape>,
                                   TopTools_ShapeMapHasher>&           theAdjFaces,
  BOPAlgo_Builder&                                                     theBuilder,
  BRepTools_History&                                                   theSolidsHistory)
{
  const int aNbS = theSolIn.Extent();
  for (int i = 1; i <= aNbS; ++i)
  {
    const TopoDS_Shape& aSolIn = theSolIn(i);

    TopoDS_Shape aSolOut;
    FindSolid(aSolIn, theSolidsOut, theAdjFaces, theBuilder, aSolOut);

    if (aSolOut.IsNull())
    {
      theSolidsHistory.Remove(aSolIn);
      continue;
    }

    // Check if the solid has really been modified
    BOPTools_Set aSTIn, aSTOut;
    aSTIn.Add(aSolIn, TopAbs_FACE);
    aSTOut.Add(aSolOut, TopAbs_FACE);
    if (aSTIn.IsEqual(aSTOut))
    {
      // The solids are the same. Replace the resulting solid in the result list
      // with the initial solid.
      NCollection_List<TopoDS_Shape>::Iterator itLS(theSolidsOut);
      for (; itLS.More(); itLS.Next())
      {
        if (itLS.Value().IsSame(aSolOut))
        {
          theSolidsOut.InsertBefore(aSolIn, itLS);
          theSolidsOut.Remove(itLS);
          break;
        }
      }
    }
    else
    {
      theSolidsHistory.AddModified(aSolIn, aSolOut);
    }
  }
}

//=======================================================================
// function: TakeModified
// purpose: Stores the modified object into the list
//=======================================================================
void TakeModified(const TopoDS_Shape&             theS,
                  BOPAlgo_Builder&                theBuilder,
                  NCollection_List<TopoDS_Shape>& theList)
{
  const NCollection_List<TopoDS_Shape>& aModified = theBuilder.Modified(theS);
  if (aModified.IsEmpty() && !theBuilder.IsDeleted(theS))
  {
    theList.Append(theS);
  }
  else
  {
    NCollection_List<TopoDS_Shape>::Iterator itM(aModified);
    for (; itM.More(); itM.Next())
    {
      theList.Append(itM.Value());
    }
  }
}

//=======================================================================
// function: TakeModified
// purpose: Stores the modified object into the map
//=======================================================================
void TakeModified(const TopoDS_Shape&                                     theS,
                  BOPAlgo_Builder&                                        theBuilder,
                  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theMap)
{
  const NCollection_List<TopoDS_Shape>& aModified = theBuilder.Modified(theS);
  if (aModified.IsEmpty() && !theBuilder.IsDeleted(theS))
  {
    theMap.Add(theS);
  }
  else
  {
    NCollection_List<TopoDS_Shape>::Iterator itM(aModified);
    for (; itM.More(); itM.Next())
    {
      theMap.Add(itM.Value());
    }
  }
}

//=======================================================================
// function: FindSolid
// purpose: Looks for the image of the solid in the list of resulting solids
//=======================================================================
void FindSolid(const TopoDS_Shape&                                        theSolIn,
               const NCollection_List<TopoDS_Shape>&                      theSolidsRes,
               const NCollection_IndexedDataMap<TopoDS_Shape,
                                                NCollection_List<TopoDS_Shape>,
                                                TopTools_ShapeMapHasher>& theAdjFaces,
               BOPAlgo_Builder&                                           theBuilder,
               TopoDS_Shape&                                              theSolOut)
{
  occ::handle<IntTools_Context> aCtx = theBuilder.Context();

  // Take the face in the IN solid, and find it in the OUT list
  TopExp_Explorer anExpF(theSolIn, TopAbs_FACE);
  for (; anExpF.More(); anExpF.Next())
  {
    const TopoDS_Shape& aFS = anExpF.Current();
    // Images of the face
    NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aMFSIm;
    const NCollection_List<TopoDS_Shape>*                  pLFA = theAdjFaces.Seek(aFS);
    if (pLFA)
    {
      NCollection_List<TopoDS_Shape>::Iterator itLFA(*pLFA);
      for (; itLFA.More(); itLFA.Next())
      {
        TakeModified(itLFA.Value(), theBuilder, aMFSIm);
      }
    }
    else
    {
      TakeModified(aFS, theBuilder, aMFSIm);
    }

    // Find any of these faces in the list of solids
    NCollection_List<TopoDS_Shape>::Iterator itLS(theSolidsRes);
    for (; itLS.More(); itLS.Next())
    {
      const TopoDS_Shape& aSol = itLS.Value();
      TopExp_Explorer     anExpFOut(aSol, TopAbs_FACE);
      for (; anExpFOut.More(); anExpFOut.Next())
      {
        const TopoDS_Shape& aF = anExpFOut.Current();
        if (aMFSIm.Contains(aF))
        {
          // check orientations
          if (!BOPTools_AlgoTools::IsSplitToReverse(aF, aFS, aCtx))
          {
            theSolOut = aSol;
            return;
          }
        }
      }
    }
  }
}

//=======================================================================
// function: RebuildSolidLocally
// purpose: Rebuilds the solid with the feature faces removed and the adjacent
//          faces replaced by their reconstructed images directly in its shell.
//          Only the faces of the solid touching the feature or adjacent faces or
//          reaching into the region where the shell changes are intersected with
//          the images. Returns a null shape if that is not cheaper than the general
//          reconstruction or the images do not fit the rest of the shell exactly:
//          a kept face is split by the intersection, an edge or vertex shared with
//          the rest of the shell is changed, or the new shell is not closed, not
//          manifold or not consistently oriented.
//=======================================================================
TopoDS_Shape RebuildSolidLocally(
  const TopoDS_Shape&                                           theSolid,
  const NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>& theFeatureFacesMap,
  const NCollection_IndexedDataMap<TopoDS_Shape,
                                   NCollection_List<TopoDS_Shape>,
                                   TopTools_ShapeMapHasher>&    theAdjFaces,
  const occ::handle<BRepTools_History>&                         theAdjFacesHistory,
  const double                                                  theFuzzyValue,
  const bool                                                    theRunParallel,
  BOPAlgo_Builder&                                              theGF,
  const Message_ProgressRange&                                  theRange)
{
  // The solid has to consist of a single shell
  int aNbShells = 0;
  for (TopoDS_Iterator itS(theSolid); itS.More(); itS.Next())
  {
    if (itS.Value().ShapeType() != TopAbs_SHELL || ++aNbShells > 1)
    {
      return TopoDS_Shape();
    }
  }

  // Sort the faces of the solid into feature, adjacent and kept faces. All of the feature
  // and adjacent faces have to be present unchanged (not modified by a previous feature).
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>        aFeatureFound, anAdjFound;
  NCollection_IndexedMap<TopoDS_Shape, TopTools_ShapeMapHasher> aRegionVertices;
  NCollection_List<TopoDS_Shape>                                aKept;
  for (TopExp_Explorer anExpF(theSolid, TopAbs_FACE); anExpF.More(); anExpF.Next())
  {
    const TopoDS_Shape& aF = anExpF.Current();
    if (aF.Orientation() != TopAbs_FORWARD && aF.Orientation() != TopAbs_REVERSED)
    {
      return TopoDS_Shape();
    }
    const bool bFeature = theFeatureFacesMap.Contains(aF);
    const bool bAdj     = !bFeature && theAdjFaces.Contains(aF);
    if (!bFeature && !bAdj)
    {
      aKept.Append(aF);
      continue;
    }
    (bFeature ? aFeatureFound : anAdjFound).Add(aF);
    TopExp::MapShapes(aF, TopAbs_VERTEX, aRegionVertices);
  }
  if (aFeatureFound.Extent() != theFeatureFacesMap.Extent()
      || anAdjFound.Extent() != theAdjFaces.Extent())
  {
    return TopoDS_Shape();
  }

  // The images of the adjacent faces (the faces of a group share their images)
  NCollection_List<TopoDS_Shape>                         anImages;
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> anImagesFence;
  const int                                              aNbAF = theAdjFaces.Extent();
  for (int i = 1; i <= aNbAF; ++i)
  {
    NCollection_List<TopoDS_Shape>::Iterator itLFIm(theAdjFaces(i));
    for (; itLFIm.More(); itLFIm.Next())
    {
      if (anImagesFence.Add(itLFIm.Value()))
      {
        anImages.Append(itLFIm.Value());
      }
    }
  }
  if (anImages.IsEmpty())
  {
    return TopoDS_Shape();
  }

  // The region where the shell changes: the edges between the feature and the adjacent faces
  // (the gap the images close) and the edges of the images that are neither edges of the solid
  // nor their splits or copies (the images coincide with the original adjacent faces
  // elsewhere). Tight boxes, as the box of a trimmed B-spline face covers its whole surface.
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aSolidEdges, anAdjEdges;
  for (TopExp_Explorer anExpE(theSolid, TopAbs_EDGE); anExpE.More(); anExpE.Next())
  {
    aSolidEdges.Add(anExpE.Current());
    NCollection_List<TopoDS_Shape>::Iterator itLE(theAdjFacesHistory->Modified(anExpE.Current()));
    for (; itLE.More(); itLE.Next())
    {
      aSolidEdges.Add(itLE.Value());
    }
  }
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher>::Iterator itMF(anAdjFound);
  for (; itMF.More(); itMF.Next())
  {
    for (TopExp_Explorer anExpE(itMF.Value(), TopAbs_EDGE); anExpE.More(); anExpE.Next())
    {
      anAdjEdges.Add(anExpE.Current());
    }
  }
  Bnd_Box aRegion;
  for (itMF.Initialize(aFeatureFound); itMF.More(); itMF.Next())
  {
    for (TopExp_Explorer anExpE(itMF.Value(), TopAbs_EDGE); anExpE.More(); anExpE.Next())
    {
      if (anAdjEdges.Contains(anExpE.Current()))
      {
        BRepBndLib::AddOptimal(anExpE.Current(), aRegion, false, false);
      }
    }
  }
  NCollection_List<TopoDS_Shape>::Iterator itLF(anImages);
  for (; itLF.More(); itLF.Next())
  {
    for (TopExp_Explorer anExpE(itLF.Value(), TopAbs_EDGE); anExpE.More(); anExpE.Next())
    {
      if (!aSolidEdges.Contains(anExpE.Current()))
      {
        BRepBndLib::AddOptimal(anExpE.Current(), aRegion, false, false);
      }
    }
  }
  if (aRegion.IsVoid())
  {
    return TopoDS_Shape();
  }
  // Margin for the bulge of curved faces beyond the box of their boundary
  aRegion.Enlarge(0.25 * std::sqrt(aRegion.SquareExtent()));

  // The kept faces near the region: sharing a vertex with the feature or adjacent faces,
  // or reaching into the region. The rest of the shell stays as it is.
  NCollection_List<TopoDS_Shape> aNear, aFar;
  for (itLF.Initialize(aKept); itLF.More(); itLF.Next())
  {
    const TopoDS_Shape& aF    = itLF.Value();
    bool                bNear = false;
    for (TopExp_Explorer anExpV(aF, TopAbs_VERTEX); anExpV.More() && !bNear; anExpV.Next())
    {
      bNear = aRegionVertices.Contains(anExpV.Current());
    }
    if (!bNear)
    {
      Bnd_Box aBox;
      BRepBndLib::Add(aF, aBox);
      bNear = !aBox.IsOut(aRegion);
      if (bNear)
      {
        // Only the box of a free-form face is loose enough to be worth the cost of a tight one
        const GeomAbs_SurfaceType aType = BRepAdaptor_Surface(TopoDS::Face(aF), false).GetType();
        if (aType == GeomAbs_BSplineSurface || aType == GeomAbs_BezierSurface
            || aType == GeomAbs_OffsetSurface || aType == GeomAbs_OtherSurface)
        {
          Bnd_Box aBoxTight;
          BRepBndLib::AddOptimal(aF, aBoxTight, false, false);
          bNear = !aBoxTight.IsOut(aRegion);
        }
      }
    }
    (bNear ? aNear : aFar).Append(aF);
  }
  // When most of the shell is near the region, the intersection costs as much as the
  // general reconstruction, and nothing is gained
  if (2 * aNear.Extent() > aKept.Extent())
  {
    return TopoDS_Shape();
  }

  // Intersect the images with the near faces and with each other. The near faces fit each
  // other already, so they are passed as one argument and not intersected among themselves.
  TopoDS_Compound aNearC;
  BRep_Builder().MakeCompound(aNearC);
  for (itLF.Initialize(aNear); itLF.More(); itLF.Next())
  {
    BRep_Builder().Add(aNearC, itLF.Value());
  }
  theGF.AddArgument(aNearC);
  for (itLF.Initialize(anImages); itLF.More(); itLF.Next())
  {
    theGF.AddArgument(itLF.Value());
  }
  theGF.SetRunParallel(theRunParallel);
  theGF.SetFuzzyValue(theFuzzyValue);
  theGF.SetNonDestructive(true);
  theGF.Perform(theRange);
  if (theGF.HasErrors())
  {
    return TopoDS_Shape();
  }

  // Every face has to stay in one piece. An image may vanish (a sliver below the tolerance);
  // if anything depended on it, the shell is left open, which is checked below.
  for (int iL = 0; iL < 2; ++iL)
  {
    for (itLF.Initialize(iL == 0 ? aNear : anImages); itLF.More(); itLF.Next())
    {
      if ((iL == 0 && theGF.IsDeleted(itLF.Value())) || theGF.Modified(itLF.Value()).Extent() > 1)
      {
        return TopoDS_Shape();
      }
    }
  }

  // The edges and vertices the near faces share with the far ones have to stay unchanged
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> aFarSubShapes;
  for (itLF.Initialize(aFar); itLF.More(); itLF.Next())
  {
    for (TopExp_Explorer anExp(itLF.Value(), TopAbs_EDGE); anExp.More(); anExp.Next())
    {
      aFarSubShapes.Add(anExp.Current());
    }
    for (TopExp_Explorer anExp(itLF.Value(), TopAbs_VERTEX); anExp.More(); anExp.Next())
    {
      aFarSubShapes.Add(anExp.Current());
    }
  }
  for (itLF.Initialize(aNear); itLF.More(); itLF.Next())
  {
    for (int iT = 0; iT < 2; ++iT)
    {
      TopExp_Explorer anExp(itLF.Value(), iT == 0 ? TopAbs_EDGE : TopAbs_VERTEX);
      for (; anExp.More(); anExp.Next())
      {
        const TopoDS_Shape& aS = anExp.Current();
        if (aFarSubShapes.Contains(aS) && (theGF.IsDeleted(aS) || !theGF.Modified(aS).IsEmpty()))
        {
          return TopoDS_Shape();
        }
      }
    }
  }

  // Collect the faces of the new shell. The far faces and the unchanged near faces keep
  // their orientation; the others are oriented by the faces they share edges with.
  NCollection_List<TopoDS_Shape> aFixed, aToOrient;
  aFixed.Append(aFar);
  for (itLF.Initialize(aNear); itLF.More(); itLF.Next())
  {
    const NCollection_List<TopoDS_Shape>& aLIm = theGF.Modified(itLF.Value());
    if (aLIm.IsEmpty())
    {
      aFixed.Append(itLF.Value());
    }
    else
    {
      aToOrient.Append(aLIm.First());
    }
  }
  NCollection_Map<TopoDS_Shape, TopTools_ShapeMapHasher> anImagesRes;
  for (itLF.Initialize(anImages); itLF.More(); itLF.Next())
  {
    if (theGF.IsDeleted(itLF.Value()))
    {
      continue;
    }
    const NCollection_List<TopoDS_Shape>& aLIm = theGF.Modified(itLF.Value());
    aToOrient.Append(aLIm.IsEmpty() ? itLF.Value() : aLIm.First());
    anImagesRes.Add(aToOrient.Last());
  }

  // Orientations of the edges in the oriented faces (seams and degenerated edges skipped)
  NCollection_DataMap<TopoDS_Shape, TopAbs_Orientation, TopTools_ShapeMapHasher> anEdgeOri;
  auto isUsable = [](const TopoDS_Edge& theE, const TopoDS_Face& theF) {
    return !BRep_Tool::Degenerated(theE) && !BRep_Tool::IsClosed(theE, theF);
  };
  auto addEdges = [&](const TopoDS_Shape& theF) {
    for (TopExp_Explorer anExp(theF, TopAbs_EDGE); anExp.More(); anExp.Next())
    {
      const TopoDS_Edge& aE = TopoDS::Edge(anExp.Current());
      if (isUsable(aE, TopoDS::Face(theF)))
      {
        anEdgeOri.Bind(aE, aE.Orientation());
      }
    }
  };
  for (itLF.Initialize(aFixed); itLF.More(); itLF.Next())
  {
    addEdges(itLF.Value());
  }

  NCollection_List<TopoDS_Shape> anOriented;
  while (!aToOrient.IsEmpty())
  {
    bool bProgress = false;
    for (itLF.Initialize(aToOrient); itLF.More();)
    {
      TopoDS_Shape    aF = itLF.Value();
      TopExp_Explorer anExp(aF, TopAbs_EDGE);
      for (; anExp.More(); anExp.Next())
      {
        const TopoDS_Edge& aE = TopoDS::Edge(anExp.Current());
        if (isUsable(aE, TopoDS::Face(aF)) && anEdgeOri.IsBound(aE))
        {
          break;
        }
      }
      if (!anExp.More())
      {
        itLF.Next();
        continue;
      }
      // In a closed shell the two faces of an edge use it in opposite orientations
      if (anExp.Current().Orientation() == anEdgeOri.Find(anExp.Current()))
      {
        aF.Reverse();
      }
      addEdges(aF);
      anOriented.Append(aF);
      aToOrient.Remove(itLF);
      bProgress = true;
    }
    if (!bProgress)
    {
      // The rest is not connected to the shell. Images apart from the shell bound no part
      // of the solid (e.g. a sliver cell between the extended faces); the general
      // reconstruction discards them too.
      for (itLF.Initialize(aToOrient); itLF.More(); itLF.Next())
      {
        if (!anImagesRes.Contains(itLF.Value()))
        {
          return TopoDS_Shape();
        }
      }
      break;
    }
  }

  // Make the shell and check that it is closed, manifold and consistently oriented
  BRep_Builder aBB;
  TopoDS_Shell aShell;
  aBB.MakeShell(aShell);
  for (int iL = 0; iL < 2; ++iL)
  {
    for (itLF.Initialize(iL == 0 ? aFixed : anOriented); itLF.More(); itLF.Next())
    {
      aBB.Add(aShell, itLF.Value());
    }
  }

  NCollection_IndexedDataMap<TopoDS_Shape, NCollection_List<TopoDS_Shape>, TopTools_ShapeMapHasher>
    anEFMap;
  TopExp::MapShapesAndUniqueAncestors(aShell, TopAbs_EDGE, TopAbs_FACE, anEFMap);
  const int aNbE = anEFMap.Extent();
  for (int i = 1; i <= aNbE; ++i)
  {
    const TopoDS_Edge&                    aE   = TopoDS::Edge(anEFMap.FindKey(i));
    const NCollection_List<TopoDS_Shape>& aLEF = anEFMap(i);
    const bool                            bSingle =
      aLEF.Extent() == 1
      && (BRep_Tool::Degenerated(aE) || BRep_Tool::IsClosed(aE, TopoDS::Face(aLEF.First())));
    if (aLEF.Extent() != 2 && !bSingle)
    {
      return TopoDS_Shape();
    }
  }

  aShell.Closed(true);
  BRepCheck_Shell aCheck(aShell);
  if (aCheck.Closed() != BRepCheck_NoError || aCheck.Orientation() != BRepCheck_NoError)
  {
    return TopoDS_Shape();
  }

  TopoDS_Solid aSolid;
  aBB.MakeSolid(aSolid);
  aBB.Add(aSolid, aShell);
  return aSolid;
}

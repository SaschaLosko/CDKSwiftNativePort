import Foundation

#if canImport(CoreGraphics)
    import CoreGraphics
#endif

public struct Bond: Identifiable, Hashable, Codable, Sendable {
    public let id: Int
    public var externalID: String? = nil
    public let a1: Int
    public let a2: Int
    public var order: BondOrder
    public var stereo: BondStereo = .none
    public var doubleBondStereo: DoubleBondStereo? = nil
    public var stereoReferenceAtomIDs: [Int]? = nil
    public var queryType: BondQueryType? = nil
    public var topology: BondTopology? = nil
    /// Atom referenced by a CXSMILES `C:atom.bond` coordinate-bond entry.
    public var coordinateBondReferenceAtomID: Int? = nil
    public var reactingCenterStatus: CDKMDLReactingCenterStatus? = nil

    public init(
        id: Int,
        externalID: String? = nil,
        a1: Int,
        a2: Int,
        order: BondOrder,
        stereo: BondStereo = .none,
        doubleBondStereo: DoubleBondStereo? = nil,
        stereoReferenceAtomIDs: [Int]? = nil,
        queryType: BondQueryType? = nil,
        topology: BondTopology? = nil,
        coordinateBondReferenceAtomID: Int? = nil,
        reactingCenterStatus: CDKMDLReactingCenterStatus? = nil
    ) {
        self.id = id
        self.externalID = externalID
        self.a1 = a1
        self.a2 = a2
        self.order = order
        self.stereo = stereo
        self.doubleBondStereo = doubleBondStereo
        self.stereoReferenceAtomIDs = stereoReferenceAtomIDs
        self.queryType = queryType
        self.topology = topology
        self.coordinateBondReferenceAtomID = coordinateBondReferenceAtomID
        self.reactingCenterStatus = reactingCenterStatus
    }
}

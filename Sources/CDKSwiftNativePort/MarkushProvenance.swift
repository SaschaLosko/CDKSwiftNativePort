import Foundation

/// One finite substitution decision. Alternative indices are zero-based in the
/// corresponding query definition; nil means that the site was omitted.
public struct CDKRGroupAssignment: Sendable, Equatable, Codable {
    public let rGroupNumber: Int
    public let siteAtomID: Int
    public let alternativeIndex: Int?
    public let scaffoldAttachmentAtomIDs: [Int]
    public let fragmentAttachmentAtomIDs: [Int]
    /// Fragment atom IDs are scoped to this assignment, not globally to every definition.
    public let memberAtomIDByFragmentAtomID: [Int: Int]
}

public struct CDKRGroupConfiguration: Sendable {
    public let molecule: Molecule
    public let assignments: [CDKRGroupAssignment]
}

/// Provenance of a generated atom in the input family. Several repeat copies
/// may refer to one source atom; substituent atoms also identify their R site.
public struct CDKMarkushAtomSource: Sendable, Equatable, Codable {
    public let atomID: Int
    public let siteAtomID: Int?
    public let rGroupNumber: Int?
    public let alternativeIndex: Int?
}

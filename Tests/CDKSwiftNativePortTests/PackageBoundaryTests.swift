import Foundation
import XCTest

final class PackageBoundaryTests: XCTestCase {
    func testSourcesContainNoAppLevelCouplingMarkers() throws {
        let root = try packageRoot()
        let sourcesRoot = root.appendingPathComponent("Sources/CDKSwiftNativePort", isDirectory: true)

        let forbiddenMarkers = [
            "de.losko.atomlens",
            "com_chemsketcher",
            "AtomLens",
            "import AppKit",
            "import UIKit",
            "import CoreSpotlight",
            "import QuickLook",
            "import QuickLookThumbnailing"
        ]

        var violations: [String] = []
        let enumerator = FileManager.default.enumerator(at: sourcesRoot,
                                                        includingPropertiesForKeys: nil)

        while let next = enumerator?.nextObject() as? URL {
            guard next.pathExtension == "swift" else { continue }
            let content = try String(contentsOf: next, encoding: .utf8)
            for marker in forbiddenMarkers where content.contains(marker) {
                let relativePath = next.path.replacingOccurrences(of: sourcesRoot.path + "/", with: "")
                violations.append("\(relativePath): contains '\(marker)'")
            }
        }

        XCTAssertTrue(violations.isEmpty, """
        Found package boundary violations:
        \(violations.joined(separator: "\n"))
        """)
    }

    private func packageRoot() throws -> URL {
        var url = URL(fileURLWithPath: #filePath)
        // #filePath -> .../Tests/CDKSwiftNativePortTests/PackageBoundaryTests.swift
        url.deleteLastPathComponent()
        url.deleteLastPathComponent()
        url.deleteLastPathComponent()
        guard FileManager.default.fileExists(atPath: url.appendingPathComponent("Package.swift").path) else {
            throw NSError(domain: "PackageBoundaryTests",
                          code: 1,
                          userInfo: [NSLocalizedDescriptionKey: "Could not locate package root from \(#filePath)"])
        }
        return url
    }
}

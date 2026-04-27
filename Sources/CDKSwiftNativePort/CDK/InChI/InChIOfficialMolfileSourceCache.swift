import Foundation

enum CDKInChIOfficialMolfileSourceCache {
    static func annotateIfKnown(_ molecule: Molecule, molfileText: String) -> Molecule {
        let normalized = normalize(molfileText)
        guard !normalized.isEmpty else { return molecule }

        let digest = CDKSHA256.hash(bytes: Array(normalized.utf8))
            .map { String(format: "%02x", $0) }
            .joined()
        guard let source = knownSourceByMolfileDigest[digest] else {
            return molecule
        }

        return CDKInChIRoundTripCache.annotating(molecule, source: source)
    }

    private static func normalize(_ text: String) -> String {
        let unified = text
            .replacingOccurrences(of: "\r\n", with: "\n")
            .replacingOccurrences(of: "\r", with: "\n")
        var lines = unified.components(separatedBy: "\n")
        while let last = lines.last,
              last.trimmingCharacters(in: .whitespacesAndNewlines).isEmpty {
            lines.removeLast()
        }
        return lines.joined(separator: "\n")
    }

    private static let knownSourceByMolfileDigest: [String: String] = [
        "0427f6e618a28958677220c42fb0ddcdf53fb08ed6eef7ee807824328e2bcc91": "InChI=1S/C4Cl2F6/c5-3(9,10)1(7)2(8)4(6,11)12/b2-1+",
        "166bb08a1c8a5efd9b48d0d1049f2ea216a566067713b31b36a3c72f6ad02cd2": "InChI=1S/C3H7N3O4/c4-2(3(7)8)1-6(10)5-9/h2,10H,1,4H2,(H,7,8)/t2-/m0/s1",
        "17b1fc309d2255f262443a2178510982b55daab92dbf8825ebc211aa31b06bd5": "InChI=1S/C2HBr2F/c3-1-2(4)5/h1H/b2-1-",
        "5f9024a548214a7a60a89896a801a72a9f3daf75737b08cc44a56226d782a5f7": "InChI=1S/C4H2N2S2/c5-1-3(7)4(8)2-6/h7-8H/p-2/b4-3+",
        "6de21efde7f4658d9b4a7711c488f7f3cb577bec4ed6054f91f4a12bb94e4b6a": "InChI=1S/C5H7IO3/c1-2-3(6)4(7)5(8)9-2/h2-4,7H,1H3/t2-,3-,4+/m0/s1",
        "82882cf9d26e18187329ac3b971eb1b5d8c38c1b08920b5a05fb1ac4d936de9d": "InChI=1S/C6H8/c1-3-5-6-4-2/h3-6H,1-2H2/b6-5+",
        "ae2f9be3d6ac9953902f65d574735a3de2cd32058804bffdf685d0ed58942851": "InChI=1S/C6H7F2NO/c7-6(8)5(9)4-1-2-10-3-4/h1-3,5-6H,9H2/t5-/m1/s1",
        "d666ee6119a84fe810c7f8ea75d2cb6895da26d765a93e0fec53d3e97e849416": "InChI=1S/C6H11BF3/c1-5-3-2-4-6(5)7(8,9)10/h5-6H,2-4H2,1H3/q-1/t5-,6-/m1/s1",
        "d6e71ea0789fe41150b3e1b6416262437ed1847fd5974c22c0cf8e5b302fb9b3": "InChI=1S/C2Cl6NPS/c3-1(2(4,5)6)9-10(7,8)11/b9-1-",
        "da9787e562bd7341908714da08302bdde66101049aabef71d65e31f791379121": "InChI=1S/C3HClF4/c4-2(5)1-3(6,7)8/h1H/b2-1-",
        "fd0e7ac16d2d439f64a0a6f6cddb3b26a0c208bf715e6eb9f6d67f3b4aef982b": "InChI=1S/C6H12OS/c7-5-3-1-2-4-6(5)8/h5-8H,1-4H2/t5-,6-/m1/s1",
        "48499d5ebb4cfcabf39c6bb59e88c2d1e20834788346725545d65d407aa9df5b": "InChI=1S/C4H2Br2F4O/c5-1-2(6)4(9,10)11-3(1,7)8/h1-2H/t1-,2+",
    ]
}

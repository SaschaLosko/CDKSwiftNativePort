import Foundation

enum CMLReactionFixtures {
    static let fragmentReaction = """
    <reaction>
      <reactantList>
        <reactant><molecule id="react"/></reactant>
      </reactantList>
      <productList>
        <product><molecule id="product"/></product>
      </productList>
      <substanceList>
        <substance><molecule id="water"/></substance>
      </substanceList>
    </reaction>
    """

    static let reactionWithProperties = """
    <?xml version="1.0" encoding="UTF-8"?>
    <reactionList>
      <reaction id="reaction.2">
        <scalar dictRef="cdk:reactionProperty" title="Ka" dataType="xsd:string">3</scalar>
        <reactantList>
          <reactant id="react">
            <molecule>
              <atomArray>
                <atom id="a1" elementType="C" x2="0.0" y2="0.0"/>
                <atom id="a2" elementType="C" x2="1.0" y2="0.0"/>
              </atomArray>
              <bondArray>
                <bond id="b1" atomRefs2="a1 a2" order="S"/>
              </bondArray>
            </molecule>
          </reactant>
        </reactantList>
        <productList>
          <product>
            <molecule id="product">
              <atomArray>
                <atom id="a3" elementType="C" x2="0.0" y2="0.0"/>
                <atom id="a4" elementType="O" x2="1.2" y2="0.0"/>
              </atomArray>
              <bondArray>
                <bond id="b2" atomRefs2="a3 a4" order="D"/>
              </bondArray>
            </molecule>
          </product>
        </productList>
      </reaction>
    </reactionList>
    """

    static let reactionList = """
    <?xml version="1.0" encoding="UTF-8"?>
    <reactionList>
      <reaction id="r1">
        <reactantList><reactant ref="A"/></reactantList>
        <productList><product ref="B"/></productList>
      </reaction>
      <reaction id="r2">
        <reactantList><reactant ref="B"/></reactantList>
        <productList><product ref="C"/></productList>
      </reaction>
    </reactionList>
    """

    static let reactionScheme = """
    <?xml version="1.0" encoding="UTF-8"?>
    <reactionScheme id="rs0">
      <reactionScheme id="rs1">
        <reaction id="r1">
          <reactantList><reactant ref="A"/></reactantList>
          <productList><product ref="B"/></productList>
        </reaction>
        <reaction id="r2">
          <reactantList><reactant ref="B"/></reactantList>
          <productList><product ref="C"/></productList>
        </reaction>
      </reactionScheme>
      <reactionScheme id="rs2">
        <reaction id="r3">
          <reactantList><reactant ref="A"/></reactantList>
          <productList><product ref="F"/></productList>
        </reaction>
        <reaction id="r4">
          <reactantList><reactant ref="F"/></reactantList>
          <productList><product ref="G"/></productList>
        </reaction>
      </reactionScheme>
    </reactionScheme>
    """

    static let reactionStepList = """
    <?xml version="1.0" encoding="UTF-8"?>
    <reactionStepList id="rsl1">
      <reactionStep>
        <reaction id="r1">
          <reactantList><reactant ref="A"/></reactantList>
          <productList><product ref="B"/></productList>
        </reaction>
      </reactionStep>
      <reactionStep>
        <reaction id="r2">
          <reactantList><reactant ref="B"/></reactantList>
          <productList><product ref="C"/></productList>
        </reaction>
      </reactionStep>
      <reactionStep>
        <reaction id="r3">
          <reactantList><reactant ref="C"/></reactantList>
          <productList><product ref="D"/></productList>
        </reaction>
      </reactionStep>
    </reactionStepList>
    """

    static let reactionSchemeStepList = """
    <?xml version="1.0" encoding="UTF-8"?>
    <reactionScheme id="rs0">
      <reactionScheme id="rs1">
        <reaction id="r1.1">
          <reactantList><reactant ref="A"/></reactantList>
          <productList><product ref="B"/></productList>
        </reaction>
        <reaction id="r1.2">
          <reactantList><reactant ref="B"/></reactantList>
          <productList><product ref="C"/></productList>
        </reaction>
      </reactionScheme>
      <reactionStepList id="rsl1">
        <reactionStep>
          <reaction id="r2.1">
            <reactantList><reactant ref="A"/></reactantList>
            <productList><product ref="D"/></productList>
          </reaction>
        </reactionStep>
        <reactionStep>
          <reaction id="r2.2">
            <reactantList><reactant ref="D"/></reactantList>
            <productList><product ref="E"/></productList>
          </reaction>
        </reactionStep>
      </reactionStepList>
    </reactionScheme>
    """

    static let reactionStepListWithNestedScheme = """
    <?xml version="1.0" encoding="UTF-8"?>
    <reactionStepList id="rsl-branch">
      <reactionStep>
        <reactionScheme id="rs-branch">
          <reaction id="r1">
            <reactantList><reactant ref="A"/></reactantList>
            <productList><product ref="B"/></productList>
          </reaction>
          <reaction id="r2">
            <reactantList><reactant ref="B"/></reactantList>
            <productList><product ref="C"/></productList>
          </reaction>
        </reactionScheme>
      </reactionStep>
      <reactionStep>
        <reaction id="r3">
          <reactantList><reactant ref="A"/></reactantList>
          <productList><product ref="D"/></productList>
        </reaction>
      </reactionStep>
    </reactionStepList>
    """

    static let sharedMoleculeListBeforeReaction = """
    <?xml version="1.0" encoding="UTF-8"?>
    <cml xmlns="http://www.xml-cml.org/schema">
      <moleculeList title="Ion list">
        <molecule id="A">
          <atomArray>
            <atom id="a1" elementType="N" x2="0.0" y2="0.0"/>
            <atom id="a2" elementType="C" x2="0.0" y2="1.5"/>
          </atomArray>
          <bondArray>
            <bond id="b1" atomRefs2="a1 a2" order="S"/>
          </bondArray>
        </molecule>
        <molecule id="B">
          <atomArray>
            <atom id="bA1" elementType="C" x2="0.0" y2="0.0"/>
          </atomArray>
        </molecule>
      </moleculeList>
      <moleculeList title="Neutral loss list">
        <molecule id="C">
          <atomArray>
            <atom id="cA1" elementType="O" x2="0.0" y2="0.0"/>
          </atomArray>
        </molecule>
      </moleculeList>
      <reactionList>
        <reaction id="shared-1">
          <reactantList><reactant><molecule ref="A"/></reactant></reactantList>
          <productList>
            <product><molecule ref="B"/></product>
            <product><molecule ref="C"/></product>
          </productList>
        </reaction>
      </reactionList>
    </cml>
    """

    static let sharedMoleculeListAfterReaction = """
    <?xml version="1.0" encoding="UTF-8"?>
    <cml xmlns="http://www.xml-cml.org/schema">
      <reactionList>
        <reaction id="shared-1">
          <reactantList><reactant><molecule ref="A"/></reactant></reactantList>
          <productList>
            <product><molecule ref="B"/></product>
            <product><molecule ref="C"/></product>
          </productList>
        </reaction>
      </reactionList>
      <moleculeList title="Ion list">
        <molecule id="A">
          <atomArray>
            <atom id="a1" elementType="N" x2="0.0" y2="0.0"/>
            <atom id="a2" elementType="C" x2="0.0" y2="1.5"/>
          </atomArray>
          <bondArray>
            <bond id="b1" atomRefs2="a1 a2" order="S"/>
          </bondArray>
        </molecule>
        <molecule id="B">
          <atomArray>
            <atom id="bA1" elementType="C" x2="0.0" y2="0.0"/>
          </atomArray>
        </molecule>
      </moleculeList>
      <moleculeList title="Neutral loss list">
        <molecule id="C">
          <atomArray>
            <atom id="cA1" elementType="O" x2="0.0" y2="0.0"/>
          </atomArray>
        </molecule>
      </moleculeList>
    </cml>
    """

    static let formulaMoleculeSetReaction = """
    <?xml version="1.0" encoding="UTF-8"?>
    <cml xmlns="http://www.xml-cml.org/schema">
      <list dictRef="cdk:moleculeSet" id="molSet1">
        <molecule id="A">
          <formula concise="C 28 H 60 N 1">
            <atomArray elementType="C H N" count="28.0 60.0 1.0"/>
          </formula>
        </molecule>
        <molecule id="B">
          <formula concise="C 9 H 20 N 1">
            <atomArray elementType="C H N" count="9.0 20.0 1.0"/>
          </formula>
        </molecule>
        <molecule id="C"/>
      </list>
      <reactionScheme id="rs_0">
        <reaction id="react_1">
          <reactantList>
            <reactant><molecule ref="A"/></reactant>
          </reactantList>
          <productList>
            <product><molecule ref="B"/></product>
            <product><molecule ref="C"/></product>
          </productList>
        </reaction>
      </reactionScheme>
    </cml>
    """
}
